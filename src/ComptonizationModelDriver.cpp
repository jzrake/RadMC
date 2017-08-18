#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include "LorentzBoost.hpp"
#include "Distributions.hpp"
#include "ComptonizationModelDriver.hpp"
#include "PathHelpers.hpp"




// ========================================================================
class SamplingScheme : public RandomVariable::SamplingScheme
{
public:
    SamplingScheme (std::function<double (double)> densityFunction, double x0, double x1)
    {
        const int numberOfTableEntries = 256;
        const double accuracyParameter = 1e-10;
        const GaussianQuadrature gauss (8);

        tabulatedCDF = TabulatedFunction::createTabulatedIntegral (
            densityFunction, x0, x1, numberOfTableEntries,
            TabulatedFunction::useEqualBinWidthsLogarithmic, gauss,
            accuracyParameter, true);
    }

    double generate (double F) override
    {
        return tabulatedCDF.lookupArgumentValue (F);
    }

private:
    TabulatedFunction tabulatedCDF;
};




// ============================================================================
ComptonizationModelDriver::ComptonizationModelDriver()
{
    addTimeSeries ("simulationTime");
    addTimeSeries ("meanPhotonEnergy");
}

void ComptonizationModelDriver::makeUserParameters (Variant::NamedValues& params)
{
    params["outdir"] = ".";
    params["tmax"] = 1.0;
    params["cpi"] = 1.0; // Checkpoint interval
    params["tsi"] = 1.0; // Time series interval
    params["theta"] = 0.01; // electron temperature in units of electron rest mass
    params["tscat"] = 1.0; // photon scattering time
    params["nphot"] = 4; // log10 of photon number
    params["nelec"] = 6; // log10 of electron number (currently for diagnostics only)
    params["ephot"] = 0.1;
    params["phot_dist_type"] = "delta"; // [delta, plaw, wein]
    params["phot_dist_lower"] = 1e-2; // in units of the mean energy
    params["phot_dist_upper"] = 1e+2; // in units of the mean energy
    params["phot_dist_index"] = -1.0; // index of the photon distribution dn/de
    params["elec_dist_type"] = "thermal"; // [thermal, plaw]
    params["elec_dist_lower"] = 1e-2; // in units of the mean energy
    params["elec_dist_upper"] = 1e+2; // in units of the mean energy
    params["elec_dist_index"] = -1.0; // index of the photon distribution dn/de
}

void ComptonizationModelDriver::configureFromParameters()
{
    // Ensure the output directory exists
    // ------------------------------------------------------------------------
    PathHelpers::ensureDirectoryExists (getParameter ("outdir"));


    // Configure the initial electron distribution
    // ------------------------------------------------------------------------
    std::string eSpectrumType = getParameter ("elec_dist_type");

    if (eSpectrumType == "thermal")
    {
        double kT = getParameter ("theta");
        double E0 = getParameter ("elec_dist_lower");
        double E1 = getParameter ("elec_dist_upper");
        double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity
        auto pdf = Distributions::makeMaxwellJuttner (kT, Distributions::Pdf);
        electronGammaBeta = RandomVariable (new SamplingScheme (pdf, E0 * urms, E1 * urms));
    }
    else if (eSpectrumType == "delta")
    {
        double kT = getParameter ("theta");
        double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity
        electronGammaBeta = RandomVariable::diracDelta (urms);
    }
    else if (eSpectrumType == "multidelta")
    {
        double kT = getParameter ("theta");
        double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity
        auto qnt = [=] (double F) { return F < 0.5 ? 0.1 * urms : std::sqrt (2) * urms; };
        electronGammaBeta = RandomVariable::fromQnt (qnt);
    }
    else if (eSpectrumType == "turbulent")
    {
        double kT = getParameter ("theta");
        double E0 = getParameter ("elec_dist_lower");
        double E1 = getParameter ("elec_dist_upper");
        double sig = std::sqrt (12. / 15 * kT);
        auto pdf = [=] (double u)
        {
           return std::sqrt (u) * std::exp (-u / sig);
        };
        electronGammaBeta = RandomVariable (new SamplingScheme (pdf, E0 * sig, E1 * sig));
    }
    else if (eSpectrumType == "plaw")
    {
        double kT = getParameter ("theta");
        double E0 = getParameter ("elec_dist_lower");
        double E1 = getParameter ("elec_dist_upper");
        double pe = getParameter ("elec_dist_index");
        electronGammaBeta = RandomVariable::powerLaw (pe, kT * E0, kT * E1);
    }
    else
    {
        throw std::runtime_error ("[ComptonizationModelDriver::configureFromParameters] "
            "unknown electron spectrum '" + eSpectrumType + "'");
    }


    // Configure the initial photon distribution
    // ------------------------------------------------------------------------
    std::string pSpectrumType = getParameter ("phot_dist_type");

    if (pSpectrumType == "delta")
    {
        double ephot = getParameter ("ephot");
        photonEnergy = RandomVariable::diracDelta (ephot);
    }
    else if (pSpectrumType == "plaw")
    {
        double ephot = getParameter ("ephot");
        double E0 = getParameter ("phot_dist_lower");
        double E1 = getParameter ("phot_dist_upper");
        double pp = getParameter ("phot_dist_index");
        photonEnergy = RandomVariable::powerLaw (pp, ephot * E0, ephot * E1);
    }
    else if (pSpectrumType == "wein")
    {
        throw std::runtime_error ("[ComptonizationModelDriver::configureFromParameters] "
            "Wein spectrum not yet implemented");
    }
    else
    {
        throw std::runtime_error ("[ComptonizationModelDriver::configureFromParameters] "
            "unknown photon spectrum '" + pSpectrumType + "'");
    }


    // Populate the initial photon distribution
    // ------------------------------------------------------------------------
    int nphot = std::pow (10, int (getParameter ("nphot")));
    nextScatteringTime = RandomVariable::exponential (getParameter ("tscat"));

    for (int n = 0; n < nphot; ++n)
    {
        Photon p = Photon::sampleIsotropic (photonEnergy);

        computeNextPhotonScatteringAndParcelVelocity (p);

        // This ensures that the right number of photons will scatter on the
        // first time step.
        // p.nextScatteringTime = -getTimestep() + nextScatteringTime.sample();
        
        photons.push_back (p);
    }


    // Write diagnostics of the initial condition
    // ------------------------------------------------------------------------
    int nelec = std::pow (10, int (getParameter ("nelec")));
    std::string filenameU = makeFilename (getParameter ("outdir"), "electron-u", ".dat");
    std::string filenameE = makeFilename (getParameter ("outdir"), "electron-e", ".dat");
    std::cout << "electron four-velocity PDF -> " << filenameU << std::endl;
    std::cout << "electron energy PDF -> " << filenameE << std::endl;
    std::ofstream outU (filenameU);
    std::ofstream outE (filenameE);
    //electronGammaBeta.outputDistribution (outU, nelec);
    electronGammaBeta.outputDistribution (outE, nelec, [] (double u) { return std::sqrt (1 + u * u) - 1; });
}

double ComptonizationModelDriver::getTimestep() const
{
    return 1.0;
}

bool ComptonizationModelDriver::shouldContinue() const
{
    Status S = getStatus();
    return S.simulationTime < double (getParameter ("tmax"));
}

void ComptonizationModelDriver::advance (double dt)
{
    Status S = getStatus();
    int numScatterings = 0;

    for (auto& p : photons)
    {
        while (p.nextScatteringTime <= S.simulationTime)
        {
            Electron e = sampleElectronForScattering (p, electronGammaBeta);
            doComptonScattering (p, e);
            p.advanceToNextScatteringTime();
            //p.nextScatteringTime = p.position.getTimeComponent() + nextScatteringTime.sample();
            computeNextPhotonScatteringAndParcelVelocity (p);

            ++numScatterings;
        }
    }
}

bool ComptonizationModelDriver::shouldRecordIterationInTimeSeries() const
{
    Status S = getStatus();
    double timeSeriesInterval = getParameter ("tsi");
    return S.simulationTime >= timeSeriesInterval * S.timeSeriesSamplesSoFar;
}

double ComptonizationModelDriver::getRecordForTimeSeries (std::string name) const
{
    if (name == "simulationTime") return getStatus().simulationTime;
    if (name == "meanPhotonEnergy") return getMeanPhotonEnergy();
    assert (false);
}

bool ComptonizationModelDriver::shouldWriteOutput() const
{
    Status S = getStatus();
    double timeBetweenOutputs = getParameter ("cpi");
    return S.simulationTime >= timeBetweenOutputs * S.outputsWrittenSoFar - 1e-12;
}

void ComptonizationModelDriver::writeOutput() const
{
    std::vector<double> energies;

    for (auto& p : photons)
    {
        energies.push_back (p.momentum.getTimeComponent());
    }

    TabulatedFunction hist = TabulatedFunction::makeHistogram (
        energies,
        256, TabulatedFunction::useEqualBinWidthsLogarithmic,
        //256, TabulatedFunction::useEqualBinWidthsLinear,
        true, true, false);

    // Write a new spectrum file
    Status S = getStatus();
    std::string filename = makeFilename (getParameter ("outdir"), "spectrum", ".dat", S.outputsWrittenSoFar);
    std::ofstream photonSpectrum (filename);
    hist.outputTable (photonSpectrum);

    // Overwrite the time-series file
    std::string timeSeriesFilename = makeFilename (getParameter ("outdir"), "time-series", ".dat");
    std::ofstream tseries (timeSeriesFilename);
    writeTimeSeriesData (tseries);
}

Electron ComptonizationModelDriver::sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta) const
{
    // Sample the electron speed (converted from gammaBeta).
    double electronU = electronGammaBeta.sample();
    double electronV = FourVector::betaFromGammaBeta (electronU);

    // Sample the electron propagation direction. The distribution of electron
    // velocities, given that scattering has occurred, is axisymmetric around
    // the photon propagation vector k, but heavier at -k than +k due to the
    // relative motion.
    RandomVariable relativeMu (Distributions::makePitchAngle (electronV, Distributions::Qnt));
    UnitVector nhat = photon.momentum.getUnitThreeVector().sampleAxisymmetric (relativeMu);

    return FourVector::fromBetaAndUnitVector (electronV, nhat);
}

void ComptonizationModelDriver::computeNextPhotonScatteringAndParcelVelocity (Photon& photon) const
{
    double betaEllStar = 0.5; // This will be taken from the turbulent field soon
    auto uf = FourVector::fromBetaAndUnitVector (betaEllStar, UnitVector::sampleIsotropic());
    double betaDotNhat = uf.getThreeVelocityAlong (photon.momentum.getUnitThreeVector());
    double effectiveTau = double (getParameter ("tscat")) * (1 - betaDotNhat) * uf.getLorentzFactor();
    photon.fluidParcelFourVelocity = uf;
    photon.nextScatteringTime = RandomVariable::exponential (effectiveTau).sample(); // exponential is cheap to create, no table
}

void ComptonizationModelDriver::doComptonScattering (Photon& photon, Electron& electron) const
{
    // Photon four-momentum in the electron rest frame
    LorentzBoost L (electron.getFourVelocity());
    FourVector p0 = photon.momentum.transformedBy (L);

    // Photon scattering angle and direction (k1) in the electron rest frame
    RandomVariable crossSection = RandomVariable::uniformOver (-1, 1); // Thomson is ~ isotropic
    double cosTheta = crossSection.sample();
    UnitVector k1 = p0.getUnitThreeVector().sampleAxisymmetric (cosTheta);

    // Photon energies in the electron frame, before and after scattering
    double e0 = p0.getTimeComponent();
    double e1 = e0 / (1 + e0 * (1 - cosTheta));

    // Photon four-momentum seen in the electron's initial rest frame (after
    // scattering) and the associated impulse - dp - seen in the lab frame.
    FourVector p1 = FourVector::nullWithUnitVector (k1) * e1;
    FourVector dp = (p1 - p0).transformedBy (L.inverted());

    // Change the photon and electron energy-momentum vector.
    photon.momentum   += dp;
    electron.momentum -= dp;
}

double ComptonizationModelDriver::getMeanPhotonEnergy() const
{
    double totalEnergy = 0.0;

    for (auto& p : photons)
    {
        totalEnergy += p.momentum.getTimeComponent();
    }
    return totalEnergy / photons.size();
}
