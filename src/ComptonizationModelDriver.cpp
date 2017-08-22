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
    addTimeSeries ("specificPhotonEnergy");
    addTimeSeries ("specificElectronEnergy");
    addTimeSeries ("specificTurbulentEnergy");
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
    params["nphot_per_mass"] = 1e-6; // Photon number per unit mass
    params["nelec"] = 6; // log10 of electron number (currently for diagnostics only)
    params["ephot"] = 0.1;
}

void ComptonizationModelDriver::configureFromParameters()
{
    // Ensure the output directory exists
    // ------------------------------------------------------------------------
    PathHelpers::ensureDirectoryExists (getParameter ("outdir"));


    // Configure the initial electron distribution
    // ------------------------------------------------------------------------
    meanParticleMass = 1836; // for Zpm = 1
    double kT = getParameter ("theta");
    regenerateElectronVelocityDistribution (kT);
    electronPopulation.momentum[0] = getSpecificInternalEnergy (kT);


    // Configure the initial photon distribution
    // ------------------------------------------------------------------------
    double ephot = getParameter ("ephot");
    int nphot = std::pow (10, int (getParameter ("nphot")));
    photonPerMass = getParameter ("nphot_per_mass"); // Photon number per unit mass
    photonEnergy = RandomVariable::diracDelta (ephot);

    for (int n = 0; n < nphot; ++n)
    {
        Photon p = Photon::sampleIsotropic (photonEnergy);
        computeNextPhotonScatteringAndParcelVelocity (p);        
        photons.push_back (p);
    }
    photonPopulation.momentum[0] = getSpecificPhotonEnergy();


    // Configure the fluid parameters
    // ------------------------------------------------------------------------
    fluidKineticEnergy = 0.0;


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
    return 0.1;
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
    double weightPerPhoton = 1.0 / photons.size();

    for (auto& p : photons)
    {
        while (p.nextScatteringTime <= S.simulationTime)
        {
            Electron e = sampleElectronForScatteringInParcel (p, electronGammaBeta);
            FourVector dp = doComptonScattering (p, e); // dp is the change in photon momentum

            photonPopulation.momentum   += dp * weightPerPhoton * photonPerMass;
            electronPopulation.momentum -= dp * weightPerPhoton * photonPerMass;

            // When dp is negative, the electron population gained energy in
            // the collision. The work done against fluid is vfluid.delta_p
            // for the photon.
            // auto vf = p.fluidParcelFourVelocity;
            // auto turbToGamma = (dp * vf + dp[0] * vf[0]) * photonPerMass;
            // electronPopulation.momentum -= dp * photonPerMass;
            // fluidKineticEnergy -= turbToGamma; // dot product of the spatial components

            ++numScatterings;

            p.advanceToNextScatteringTime();
            computeNextPhotonScatteringAndParcelVelocity (p);
        }

        // std::cout << std::setprecision(10) << "MC: " << getSpecificPhotonEnergy() << " " << photonPopulation.momentum[0] << "\n";
    }


    // Recompute the electron temperature
    // regenerateElectronVelocityDistribution (kT);

    double kT_elec = getElectronTemperature();
    double kT_phot = getMeanPhotonEnergy();

    double Eint = getSpecificInternalEnergy (kT_elec);
    double Erad = getSpecificPhotonEnergy();
    double Ekin = fluidKineticEnergy;

    std::ostringstream message;
    message << std::scientific << std::setprecision(4);
    message << "kT_elec=" << kT_elec << " ";
    message << "kT_phot=" << kT_phot << " ";
    message << "Eint=" << Eint << " ";
    message << "Erad=" << Erad << " ";
    message << "Ekin=" << Ekin << " ";
    message << "total=" << Eint + Erad + Ekin << " ";
    std::cout << message.str() << std::endl;
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
    if (name == "specificPhotonEnergy") return getSpecificPhotonEnergy();
    if (name == "specificElectronEnergy") return electronPopulation.momentum[0];
    if (name == "specificTurbulentEnergy") return fluidKineticEnergy;
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

Electron ComptonizationModelDriver::sampleElectronForScatteringInParcel (const Photon& photon, RandomVariable& electronGammaBeta) const
{
    LorentzBoost L (photon.fluidParcelFourVelocity);

    // Get the photon in the fluid rest frame
    Photon q = photon.transformedBy (L);
    Electron f = sampleElectronForScattering (q, electronGammaBeta);
    Electron e = f.transformedBy (L.inverted());
    return e;
}

void ComptonizationModelDriver::computeNextPhotonScatteringAndParcelVelocity (Photon& photon) const
{
    double betaEllStar = std::sqrt (fluidKineticEnergy); // use non-relativistic expression
    auto uf = FourVector::fromBetaAndUnitVector (betaEllStar, UnitVector::sampleIsotropic());
    double betaDotNhat = uf.getThreeVelocityAlong (photon.momentum.getUnitThreeVector());
    double effectiveTau = double (getParameter ("tscat")) * (1 - betaDotNhat) * uf.getLorentzFactor();
    photon.fluidParcelFourVelocity = uf;
    photon.nextScatteringTime = (0.0
        + photon.position.getTimeComponent()
        + RandomVariable::exponential (effectiveTau).sample()); // exponential is cheap to create, no table
}

FourVector ComptonizationModelDriver::doComptonScattering (Photon& photon, Electron& electron) const
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

    return dp;
}

double ComptonizationModelDriver::getTotalPhotonEnergyMC() const
{
    double totalEnergy = 0.0;

    for (const auto& p : photons)
    {
        totalEnergy += p.momentum[0];
    }
    return totalEnergy;
}

double ComptonizationModelDriver::getMeanPhotonEnergy() const
{
    double meanEnergyPerPhoton = getTotalPhotonEnergyMC() / photons.size();
    return meanEnergyPerPhoton;
}

double ComptonizationModelDriver::getSpecificPhotonEnergy() const
{
    double meanEnergyPerPhoton = getMeanPhotonEnergy();
    double photonEnergyPerMass = meanEnergyPerPhoton * photonPerMass;
    return photonEnergyPerMass;
}

double ComptonizationModelDriver::getElectronTemperature() const
{
    // This is actually the specific internal energy over all plasma species
    double e = electronPopulation.momentum[0];
    return 1. / 3 * meanParticleMass * e;
}

double ComptonizationModelDriver::getSpecificInternalEnergy (double electronTemperature) const
{
    return 3 * electronTemperature / meanParticleMass;
}

void ComptonizationModelDriver::regenerateElectronVelocityDistribution (double electronTemperature)
{
    double kT = electronTemperature;
    double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity

    if (kT < 1e-2)
    {
        auto pdf = Distributions::makeMaxwellBoltzmann (kT, Distributions::Pdf);
        electronGammaBeta = RandomVariable (new SamplingScheme (pdf, urms * 0.01, urms * 10));        
    }
    else
    {
        // For relativistic temperatures the sampler can tolerate larger maximum speeds
        auto pdf = Distributions::makeMaxwellJuttner (kT, Distributions::Pdf);
        electronGammaBeta = RandomVariable (new SamplingScheme (pdf, urms * 0.01, urms * 100));
    }
}
