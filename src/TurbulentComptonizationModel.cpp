#include <iostream>
#include "TurbulentComptonizationModel.hpp"
#include "Distributions.hpp"
#include "LorentzBoost.hpp"

using Photon = TurbulentComptonizationModel::Photon;
using Electron = TurbulentComptonizationModel::Electron;




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




// ========================================================================
TurbulentComptonizationModel::TurbulentComptonizationModel (Config config) : config (config)
{
    simulationTime = 0.0;

    // Configure the initial electron distribution
    // ------------------------------------------------------------------------
    double kT = config.theta;

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
    double betaTurb = double (config.beta_turb);
    fluidKineticEnergy = std::pow (betaTurb, 2);
    zPlusMinus = 1.0;


    // Configure the initial photon distribution
    // ------------------------------------------------------------------------
    double A = 1.0 / (1 + zPlusMinus / protonToElectronMassRatio);
    double ephot = config.ephot;
    int nphot = std::pow (10, int (config.nphot));
    photonMeanFreePath = double (config.ell_star);
    photonPerMass = config.nphot_per_mass; // Photon number per unit mass
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
}

double TurbulentComptonizationModel::getTimestep() const
{
    double cascadeDt = cascadeModel.getShortestTimeScale() * 0.1;
    double scatterDt = photonMeanFreePath * 0.1;
    double coolingDt = getComptonCoolingTime() * 0.1;
    return std::min ({cascadeDt, scatterDt, coolingDt});
}

void TurbulentComptonizationModel::advance (double dt)
{
    int numScatterings = 0;
    double weightPerPhoton = 1.0 / photons.size();
    double epsEllStar = cascadeModel.getEnergyFluxThroughPhotonMeanFreePathScale();

    for (auto& p : photons)
    {
        while (p.nextScatteringTime <= simulationTime)
        {
            Electron e = sampleElectronForScatteringInParcel (p, electronGammaBeta);
            FourVector dp = doComptonScattering (p, e); // dp is the change in photon momentum


            // When dp is negative, the electron population gained energy in
            // the collision. The work done against fluid is vfluid.dp.
            auto vf = p.fluidParcelFourVelocity;
            auto turbToGamma = (dp[1] * vf[1] + dp[2] * vf[2] + dp[3] * vf[3]) * weightPerPhoton * photonPerMass;
            auto turbToInt = epsEllStar * dt * weightPerPhoton; // FIX THIS!!
            auto gammaToInt = -dp[0] * weightPerPhoton * photonPerMass + turbToGamma;


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

    cascadeModel.radiativeEnergyDensity = getSpecificPhotonEnergy();
    cascadeModel.advance (dt);
    simulationTime += dt;
}

std::vector<double> TurbulentComptonizationModel::getCascadeWaveNumbers() const
{
    return cascadeModel.powerSpectrum.getDataX();
}

std::vector<double> TurbulentComptonizationModel::getCascadePowerSpectrum() const
{
    return cascadeModel.powerSpectrum.getDataY();
}

std::vector<Photon> TurbulentComptonizationModel::getPhotons() const
{
    return photons;
}

double TurbulentComptonizationModel::getFluidVelocityAtPhotonMeanFreePathScale() const
{
    // double betaEllStar = std::sqrt (fluidKineticEnergy); // use non-relativistic expression
    return cascadeModel.getEddyVelocityAtScale (photonMeanFreePath);
}

double TurbulentComptonizationModel::getTotalPhotonEnergyMC() const
{
    double totalEnergy = 0.0;

    for (const auto& p : photons)
    {
        totalEnergy += p.momentum[0];
    }
    return totalEnergy;
}

double TurbulentComptonizationModel::getMeanPhotonEnergy() const
{
    double meanEnergyPerPhoton = getTotalPhotonEnergyMC() / photons.size();
    return meanEnergyPerPhoton;
}

double TurbulentComptonizationModel::getSpecificPhotonEnergy() const
{
    double meanEnergyPerPhoton = getMeanPhotonEnergy();
    double photonEnergyPerMass = meanEnergyPerPhoton * photonPerMass;
    return photonEnergyPerMass;
}

double TurbulentComptonizationModel::getElectronTemperature() const
{
    double e = plasmaInternalEnergy;
    return 1. / 3 * meanParticleMass * e;
}

double TurbulentComptonizationModel::getPhotonTemperature() const
{
    return 1. / 3 * getMeanPhotonEnergy();
}

double TurbulentComptonizationModel::getSpecificInternalEnergy (double electronTemperature) const
{
    return 3 * electronTemperature / meanParticleMass;
}

double TurbulentComptonizationModel::getComptonCoolingTime() const
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

Electron TurbulentComptonizationModel::sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta) const
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

Electron TurbulentComptonizationModel::sampleElectronForScatteringInParcel (const Photon& photon, RandomVariable& electronGammaBeta) const
{
    LorentzBoost L (photon.fluidParcelFourVelocity);

    // Get the photon in the fluid rest frame
    Photon q = photon.transformedBy(L);
    Electron f = sampleElectronForScattering (q, electronGammaBeta);
    Electron e = f.transformedBy (L.inverted());
    return e;
}

void TurbulentComptonizationModel::computeNextPhotonScatteringAndParcelVelocity (Photon& photon) const
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

FourVector TurbulentComptonizationModel::doComptonScattering (Photon& photon, Electron& electron) const
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

void TurbulentComptonizationModel::regenerateElectronVelocityDistribution (double electronTemperature)
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
