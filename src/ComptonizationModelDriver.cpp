#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "LorentzBoost.hpp"
#include "Distributions.hpp"
#include "ComptonizationModelDriver.hpp"




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

    auto electronPdf = Distributions::makeMaxwellJuttner (kT, Distributions::Pdf);
    electronGammaBeta = RandomVariable (new SamplingScheme (electronPdf, 1e-8, 1.0));
    photonEnergy = RandomVariable::diracDelta (ephot);

    for (int n = 0; n < nphot; ++n)
    {
        photons.push_back (Photon::sampleIsotropic (photonEnergy));
    }
}

void ComptonizationModelDriver::printStartupMessage() const
{

}

double ComptonizationModelDriver::getTimestep() const
{
    return 1.0;
}

bool ComptonizationModelDriver::shouldContinue() const
{
    SimulationDriver::Status S = getStatus();
    return S.simulationTime < double (getParameter ("tmax"));
}

void ComptonizationModelDriver::advance (double dt)
{
    for (auto& p : photons)
    {
        Electron e = sampleElectronForScattering (p, electronGammaBeta);
        doComptonScattering (p, e);
    }
}

bool ComptonizationModelDriver::shouldRecordIterationInTimeSeries() const
{
    SimulationDriver::Status S = getStatus();
    double timeSeriesInterval = getParameter ("tsi");
    return S.simulationTime >= timeSeriesInterval * S.timeSeriesSamplesSoFar;
}

double ComptonizationModelDriver::getRecordForTimeSeries (std::string name) const
{
    if (name == "simulationTime") return getStatus().simulationTime;
    if (name == "meanPhotonEnergy") return getMeanPhotonEnergy();
    return 0.0;
}

bool ComptonizationModelDriver::shouldWriteOutput() const
{
    SimulationDriver::Status S = getStatus();
    double timeBetweenOutputs = getParameter ("cpi");
    return S.simulationTime >= timeBetweenOutputs * S.outputsWrittenSoFar - 1e-12;
}

void ComptonizationModelDriver::writeOutput (std::string filename) const
{
    std::vector<double> energies;

    for (auto& p : photons)
    {
        energies.push_back (p.momentum.getTimeComponent());
    }

    TabulatedFunction hist = TabulatedFunction::makeHistogram (energies, 256,
        TabulatedFunction::useEqualBinWidthsLinear, true, true, false);

    std::ofstream photonSpectrum (filename);
    hist.outputTable (photonSpectrum);
}

std::string ComptonizationModelDriver::makeOutputFilename() const
{
    SimulationDriver::Status S = getStatus();
    std::ostringstream filenameStream;
    filenameStream << getParameter ("outdir") << "/photspec.";
    filenameStream << std::setfill ('0') << std::setw (6) << S.outputsWrittenSoFar;
    filenameStream << ".dat";
    return filenameStream.str();
}

Electron ComptonizationModelDriver::sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta)
{
    // Sample the electron speed (Note: convert from gammaBeta if MJ).
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
