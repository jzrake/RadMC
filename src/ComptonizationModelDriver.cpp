#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "LorentzBoost.hpp"
#include "Distributions.hpp"
#include "ComptonizationModelDriver.hpp"
#include "PathHelpers.hpp"





namespace Cow
{
    class Timer;
}




#include <ctime>
#include <iomanip>
#include <sstream>

using namespace Cow;


/**
A handy class for checking execution time using RAII.

    auto timer = Timer();
    someExpensiveFunction();
    std::cout << timer.age() << std::endl;
*/
class Cow::Timer
{
public:
    Timer();

    /**
    Return the time, in seconds, since the timer was instantiated.
    */
    double age() const;

    /**
    Return the timer age, in minutes.
    */
    double minutes() const;

    /**
    Return the age, formatted as e.g. "21.4 seconds".
    */
    std::string ageInSeconds() const;

private:
    std::clock_t timeInstantiated;
};



// ============================================================================
Timer::Timer() : timeInstantiated (std::clock())
{

}

double Timer::age() const
{
    return double (std::clock() - timeInstantiated) / CLOCKS_PER_SEC;
}

double Timer::minutes() const
{
    return age() / 60;
}

std::string Timer::ageInSeconds() const
{
    auto stream = std::ostringstream();
    stream << std::setprecision(2) << age() << " seconds";
    return stream.str();
}

static Timer timer;



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
    addTimeSeries ("photonTemperature");
    addTimeSeries ("electronTemperature");
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
    params["ell_star"] = 1e-3; // photon mean free path
    params["nphot"] = 4; // log10 of photon number
    params["nphot_per_mass"] = 1e-1; // Photon number per unit mass (1 means mp / me photons per baryon)
    params["nelec"] = 5; // log10 of electron number (currently for diagnostics only)
    params["ephot"] = 1e-3;
    params["beta_turb"] = 0.5;
}

void ComptonizationModelDriver::configureFromParameters()
{
    // Ensure the output directory exists
    // ------------------------------------------------------------------------
    PathHelpers::ensureDirectoryExists (getParameter ("outdir"));


    // Configure the initial electron distribution
    // ------------------------------------------------------------------------
    double kT = getParameter ("theta");

    if (kT == 0)
    {
        std::cout << "Using cold electron approximation\n";
        coldElectrons = true;
    }
    else
    {
        coldElectrons = false;
    }
    protonToElectronMassRatio = 1836;
    meanParticleMass = (protonToElectronMassRatio + zPlusMinus) / (1 + zPlusMinus);
    regenerateElectronVelocityDistribution (kT);
    plasmaInternalEnergy = getSpecificInternalEnergy (kT);


    // Configure the fluid parameters
    // ------------------------------------------------------------------------
    double betaTurb = double (getParameter ("beta_turb"));
    fluidKineticEnergy = std::pow (betaTurb, 2);
    zPlusMinus = 1.0;


    // Configure the initial photon distribution
    // ------------------------------------------------------------------------
    double A = 1.0 / (1 + zPlusMinus / protonToElectronMassRatio);
    double ephot = getParameter ("ephot");
    int nphot = std::pow (10, int (getParameter ("nphot")));
    photonMeanFreePath = double (getParameter ("ell_star"));
    photonPerMass = getParameter ("nphot_per_mass"); // Photon number per unit mass
    photonEnergy = RandomVariable::diracDelta (ephot);
    photonPerProton = photonPerMass / A * protonToElectronMassRatio;

    for (int n = 0; n < nphot; ++n)
    {
        Photon p = Photon::sampleIsotropic (photonEnergy);
        computeNextPhotonScatteringAndParcelVelocity (p);
        photons.push_back (p);
    }

    if (photonMeanFreePath >= 1.0)
    {
        throw std::runtime_error ("ell_star (photon mean free path) must be < 1");
    }
    cascadeModel = RichardsonCascade (10 / photonMeanFreePath, 128);
    cascadeModel.photonMeanFreePath = photonMeanFreePath;
    cascadeModel.radiativeEnergyDensity = getSpecificPhotonEnergy();
    cascadeModel.cascadePower = std::pow (betaTurb, 3);


    double Re = cascadeModel.getReynoldsNumber (betaTurb);
    std::cout << "Optical depth is " << 1.0 / photonMeanFreePath << std::endl;
    std::cout << "Reynolds number is " << Re << std::endl;
    std::cout << "Viscous cutoff expected at " << std::pow (Re, 3. / 4) << std::endl;
    std::cout << "If optical depth tau > Re^(3/4) (=" << std::pow (Re, 3. / 4) << ")" << std::endl;
    std::cout << "Photons per proton " << photonPerProton << std::endl;


    // Write diagnostics of the initial condition
    // ------------------------------------------------------------------------
    // int nelec = std::pow (10, int (getParameter ("nelec")));
    // std::string filenameU = makeFilename (getParameter ("outdir"), "electron-u", ".dat");
    // std::string filenameE = makeFilename (getParameter ("outdir"), "electron-e", ".dat");
    // std::cout << "electron four-velocity PDF -> " << filenameU << std::endl;
    // std::cout << "electron energy PDF -> " << filenameE << std::endl;
    // std::ofstream outU (filenameU);
    // std::ofstream outE (filenameE);
    // electronGammaBeta.outputDistribution (outE, nelec, [] (double u) { return std::sqrt (1 + u * u) - 1; });
}

double ComptonizationModelDriver::getTimestep() const
{
    static int counter = 0;
    static double last_dt = 0.0;

    if (counter == 0) // so we don't compute the timeset too often
    {
        double cascadeDt = cascadeModel.getShortestTimeScale() * 0.1;
        double scatterDt = photonMeanFreePath * 0.1;
        double coolingDt = getComptonCoolingTime() * 0.1;
        last_dt = std::min ({cascadeDt, scatterDt, coolingDt});
    }
    counter = (counter + 1) % 1000;

    return last_dt;
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

    // double fluidKineticEnergyAtStart = fluidKineticEnergy;
    double weightPerPhoton = 1.0 / photons.size();
    double meanScatteringTime = photonMeanFreePath;

    for (auto& p : photons)
    {
        while (p.nextScatteringTime <= S.simulationTime)
        {
            Electron e = sampleElectronForScatteringInParcel (p, electronGammaBeta);
            FourVector dp = doComptonScattering (p, e); // dp is the change in photon momentum


            // When dp is negative, the electron population gained energy in
            // the collision. The work done against fluid is vfluid.dp.
            auto vf = p.fluidParcelFourVelocity;
            auto turbToGamma = (dp[1] * vf[1] + dp[2] * vf[2] + dp[3] * vf[3]) * weightPerPhoton * photonPerMass;
            auto turbToInt = cascadeModel.getEnergyFluxThroughScale (photonMeanFreePath) * meanScatteringTime * weightPerPhoton;
            auto gammaToInt = -dp[0] * weightPerPhoton * photonPerMass + turbToGamma;

            // std::cout << dp[1] << std::endl;

            plasmaInternalEnergy += (turbToInt + gammaToInt);
            fluidKineticEnergy   -= (turbToInt + turbToGamma);

            ++numScatterings;

            p.advanceToNextScatteringTime();
            computeNextPhotonScatteringAndParcelVelocity (p);
        }
    }



    // Recompute the electron temperature
    double kT = getElectronTemperature();

    if (kT > 0.0)
    {
        regenerateElectronVelocityDistribution (getElectronTemperature());
    }
    else if (! coldElectrons)
    {
        std::cout << "WARNING: negative electron temperature " << kT << std::endl;
    }


    // double dissipatedPower = -(fluidKineticEnergy - fluidKineticEnergyAtStart) / dt;
    // double viscosity = dissipatedPower / cascadeModel.getDissipationRatePerViscosity();

    // std::cout << "dissipatedPower: " << dissipatedPower << std::endl;
    // std::cout << "dissipationRatePerViscosity: " << cascadeModel.getDissipationRatePerViscosity() << std::endl;
    // std::cout << "viscosity needed: " << viscosity << std::endl;
    // std::cout << "photon viscosity: " << photonMeanFreePath * getSpecificPhotonEnergy() << std::endl;

    cascadeModel.radiativeEnergyDensity = getSpecificPhotonEnergy();
    cascadeModel.advance (dt);
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
    if (name == "photonTemperature") return getPhotonTemperature();
    if (name == "electronTemperature") return getElectronTemperature();
    if (name == "specificPhotonEnergy") return getSpecificPhotonEnergy();
    if (name == "specificElectronEnergy") return plasmaInternalEnergy;
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
    // std::cout << timer.ageInSeconds() << std::endl;
    // timer = Timer();

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



    // Write the cascde model data
    // ------------------------------------------------------------------------
    std::vector<std::vector<double>> columns;
    columns.push_back (cascadeModel.powerSpectrum.getDataX());
    columns.push_back (cascadeModel.powerSpectrum.getDataY());
    // columns.push_back (cascadeModel.getEddyTurnoverTime());
    // columns.push_back (cascadeModel.getDampingTime());

    std::string turbFilename = makeFilename (getParameter ("outdir"), "cascade", ".dat", S.outputsWrittenSoFar);
    std::ofstream stream (turbFilename);
    writeAsciiTable (columns, stream);
}

std::string ComptonizationModelDriver::getStatusMessage() const
{
    double kT_elec = getElectronTemperature();
    double kT_phot = getPhotonTemperature();

    double Eint = getSpecificInternalEnergy (kT_elec);
    double Erad = getSpecificPhotonEnergy();
    double Ekin = fluidKineticEnergy;

    std::ostringstream message;
    message << std::scientific << std::setprecision(1);
    message << "tc=" << getComptonCoolingTime() << " ";
    message << "Te=" << kT_elec << " ";
    message << "Tr=" << kT_phot << " ";
    message << "Ei=" << Eint << " ";
    message << "Er=" << Erad << " ";
    message << "Ek=" << Ekin << " ";
    message << "vs=" << cascadeModel.getEddyVelocityAtScale (photonMeanFreePath) << " ";
    message << "ep=" << cascadeModel.getEnergyFluxThroughScale (photonMeanFreePath) << " ";
    message << "total=" << Eint + Erad + Ekin << " ";
    return message.str();
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
    Photon q = photon.transformedBy(L);
    Electron f = sampleElectronForScattering (q, electronGammaBeta);
    Electron e = f.transformedBy (L.inverted());
    return e;
}

void ComptonizationModelDriver::computeNextPhotonScatteringAndParcelVelocity (Photon& photon) const
{
    double betaEllStar = getFluidVelocityAtPhotonMeanFreePathScale();
    auto uf = FourVector::fromBetaAndUnitVector (betaEllStar, UnitVector::sampleIsotropic());
    double betaDotNhat = uf.getThreeVelocityAlong (photon.momentum.getUnitThreeVector());

    double relativeMotionFactor = (1 - betaDotNhat) * uf.getLorentzFactor();
    double photonTravelDistance = photonMeanFreePath / relativeMotionFactor;
    double actualTravel = photonTravelDistance * RandomVariable::exponential(1).sample();

    photon.fluidParcelFourVelocity = uf;
    photon.nextScatteringTime = (photon.position.getTimeComponent() + actualTravel);
}

FourVector ComptonizationModelDriver::doComptonScattering (Photon& photon, Electron& electron) const
{
    // Photon four-momentum in the electron rest frame
    LorentzBoost L (electron.getFourVelocity());
    FourVector p0 = photon.momentum.transformedBy(L);

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

double ComptonizationModelDriver::getFluidVelocityAtPhotonMeanFreePathScale() const
{
    // double betaEllStar = std::sqrt (fluidKineticEnergy); // use non-relativistic expression
    return cascadeModel.getEddyVelocityAtScale (photonMeanFreePath);
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
    double e = plasmaInternalEnergy;
    return 1. / 3 * meanParticleMass * e;
}

double ComptonizationModelDriver::getPhotonTemperature() const
{
    return 1. / 3 * getMeanPhotonEnergy();
}

double ComptonizationModelDriver::getSpecificInternalEnergy (double electronTemperature) const
{
    return 3 * electronTemperature / meanParticleMass;
}

double ComptonizationModelDriver::getComptonCoolingTime() const
{
    if (coldElectrons)
    {
        return 1.0; // Just an arbitrary "large" number.
    }
    double kT = getElectronTemperature();
    double u2 = kT < 1 ? 3 * kT : 12 * kT * kT;
    double e = plasmaInternalEnergy * protonToElectronMassRatio;
    double eps = getMeanPhotonEnergy();
    double tscat = photonMeanFreePath;
    double tcomp = tscat * 3. / 4 / photonPerProton * zPlusMinus * e / eps / u2;
    return tcomp;
}

void ComptonizationModelDriver::regenerateElectronVelocityDistribution (double electronTemperature)
{
    if (coldElectrons)
    {
        electronGammaBeta = RandomVariable::diracDelta(0.0);
        return;
    }

    double kT = electronTemperature;
    double urms = std::sqrt (kT < 1 ? 3 * kT : 12 * kT * kT); // approximate RMS four-velocity

    if (kT < 1e-2)
    {
        auto pdf = Distributions::makeMaxwellBoltzmann (kT, Distributions::Pdf);
        electronGammaBeta = RandomVariable (new SamplingScheme (pdf, urms * 0.01, urms * 10));        
    }
    else
    {
        auto pdf = Distributions::makeMaxwellJuttner (kT, Distributions::Pdf);
        electronGammaBeta = RandomVariable (new SamplingScheme (pdf, urms * 0.01, urms * 10));
    }
}
