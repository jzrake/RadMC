#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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
        const double accuracyParameter = 1e-12;
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
    params["ephot"] = 0.1;
    params["nphot"] = 4; // log10 of photon number
}

void ComptonizationModelDriver::configureFromParameters()
{
    int nphot = std::pow (10, int (getParameter ("nphot")));
    double kT = getParameter ("theta");
    double ephot = getParameter ("ephot");
    double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity

    auto electronPdf = Distributions::makeMaxwellJuttner (kT, Distributions::Pdf);
    electronGammaBeta = RandomVariable (new SamplingScheme (electronPdf, urms / 100, urms * 5));
    photonEnergy = RandomVariable::diracDelta (ephot);
    //photonEnergy = RandomVariable::powerLaw (-2, ephot * 0.01, ephot * 100.0);
    nextScatteringTime = RandomVariable::exponential (1.0);

    for (int n = 0; n < nphot; ++n)
    {
        photons.push_back (Photon::sampleIsotropic (photonEnergy));
    }

    PathHelpers::ensureDirectoryExists (getParameter ("outdir"));
}

void ComptonizationModelDriver::printStartupMessage() const
{
    std::string filenameU = makeFilename (getParameter ("outdir"), "electron-u", ".dat");
    std::string filenameE = makeFilename (getParameter ("outdir"), "electron-e", ".dat");
    std::cout << "electron four-velocity PDF -> " << filenameU << std::endl;
    std::cout << "electron four-velocity PDF -> " << filenameE << std::endl;
    std::ofstream outU (filenameU);
    std::ofstream outE (filenameE);
    electronGammaBeta.outputDistribution (outU, 1e5);
    electronGammaBeta.outputDistribution (outE, 1e5, [] (double u) { return std::sqrt (1 + u * u) - 1; });
}

double ComptonizationModelDriver::getTimestep() const
{
    return 1.;
}

bool ComptonizationModelDriver::shouldContinue() const
{
    Status S = getStatus();
    return S.simulationTime < double (getParameter ("tmax"));
}

void ComptonizationModelDriver::advance (double dt)
{
    Status S = getStatus();

    for (auto& p : photons)
    {
        //int numScatterings = 0;

        while (p.nextScatteringTime <= S.simulationTime)
        {
            Electron e = sampleElectronForScattering (p, electronGammaBeta);
            doComptonScattering (p, e);
            p.advanceToNextScatteringTime();
            p.nextScatteringTime = p.position[0] + nextScatteringTime.sample();
            //++numScatterings;
        }

        //std::cout << numScatterings << std::endl;

        //Electron e = sampleElectronForScattering (p, electronGammaBeta);
        //doComptonScattering (p, e);
        //p.advancePosition (dt);
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
    throw std::runtime_error ("unknown time series name '" + name + "'");
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

Electron ComptonizationModelDriver::sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta)
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

Electron ComptonizationModelDriver::sampleElectronForScattering (const Photon& photon, const FourVector::Field& field)
{
    return field (photon.position);
}

void ComptonizationModelDriver::doComptonScattering (Photon& photon, Electron& electron)
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
