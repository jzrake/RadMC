#pragma once

#include "TabulatedFunction.hpp"
#include "FourVector.hpp"
#include "RichardsonCascade.hpp"




class TurbulentComptonizationModel
{
public:



    struct Config
    {
        double theta = 0.01; // electron temperature in units of electron rest mass
        double ell_star = 1e-3; // photon mean free path
        double nphot = 4; // log10 of photon number
        double nphot_per_mass = 1e-1; // Photon number per unit mass (1 means mp / me photons per baryon)
        double nelec = 5; // log10 of electron number (currently for diagnostics only)
        double ephot = 1e-3;
        double beta_turb = 0.5;
    };



    struct IterationReport
    {
        double timeAfterIteration;
        double meanScatteringsPerPhoton;
        double meanScatteringAngleWithBulk;
        double meanScatteringAngleInParcel;
        double viscousPowerBasedOnPhotons;
        double viscousPowerBasedOnCascade;
    };




    class Electron
    {
    public:
        Electron () : momentum (1, 0, 0, 0) {}
        Electron (FourVector momentum) : momentum (momentum) {}

        FourVector getFourVelocity() const
        {
            return momentum;
        }

        Electron transformedBy (const LorentzBoost& L) const
        {
            return Electron (momentum.transformedBy (L));
        }

        FourVector momentum;
    };




    class Photon
    {
    public:
        static Photon sampleIsotropic (RandomVariable& photonEnergy)
        {
            double E = photonEnergy.sample();
            return FourVector::nullWithUnitVector (UnitVector::sampleIsotropic()) * E;
        }
        Photon() : nextScatteringTime (0) {}
        Photon (FourVector momentum) : momentum (momentum), nextScatteringTime (0) {}

        Photon transformedBy (const LorentzBoost& L) const
        {
            Photon prime;
            prime.position = position.transformedBy (L);
            prime.momentum = momentum.transformedBy (L);
            prime.fluidParcelFourVelocity = fluidParcelFourVelocity.transformedBy (L);
            return prime;
        }

        /**
        Update the position four-vector based on the momentum and time step dt.
        */
        void advancePosition (double dt)
        {
            position += getDisplacement (dt);
        }

        /**
        This is just a shortcut for advancing the position to the next scattering
        time.
        */
        void advanceToNextScatteringTime()
        {
            advancePosition (nextScatteringTime - position[0]);
        }

        /**
        Compute the photon's displacement vector for a time dt.
        */
        FourVector getDisplacement (double dt)
        {
            UnitVector nhat = momentum.getUnitThreeVector();
            return FourVector (dt, nhat.getX() * dt, nhat.getY() * dt, nhat.getZ() * dt);
        }

        FourVector position;
        FourVector momentum;
        FourVector fluidParcelFourVelocity;
        double nextScatteringTime;
    };




    TurbulentComptonizationModel (Config config);

    // Driver functions
    // ========================================================================
    double getTimestep() const;
    IterationReport advance (double dt);

    // Functions that query array data
    // ========================================================================
    std::vector<Photon> getPhotons() const;
    std::vector<double> getCascadeWaveNumberBins() const;
    std::vector<double> getCascadePowerSpectrum() const;
    std::vector<double> getPhotonEnergyBins() const;
    std::vector<double> getPhotonSpectrum() const;

    // Functions that query diagnostic data
    // ========================================================================
    double getElectronTemperature() const;
    double getPhotonTemperature() const;
    double getEffectiveWaveTemperature() const;
    double getComptonCoolingTime() const;
    double getSpecificInternalEnergy() const;
    double getSpecificKineticEnergy() const;
    double getSpecificPhotonEnergy() const;
    double getEddyVelocityAtScale (double ell) const;
    
private:
    double getFluidVelocityAtPhotonMeanFreePathScale() const;
    double getTotalPhotonEnergyMC() const;
    double getMeanPhotonEnergy() const;
    double getSpecificInternalEnergyForTemp (double electronTemperature) const;
    Electron sampleElectronForScattering (const Photon& photon) const;
    Electron sampleElectronForScatteringInParcel (const Photon& photon) const;
    FourVector doComptonScattering (Photon& photon, Electron& electron) const;
    void computeNextPhotonScatteringAndParcelVelocity (Photon& photon) const;
    void regenerateElectronVelocityDistribution (double electronTemperature);
    TabulatedFunction computePhotonSpectrum() const;

    Config config;
    RandomVariable electronGammaBeta;
    RandomVariable photonEnergy;
    RichardsonCascade cascadeModel;
    bool coldElectrons;
    double simulationTime;
    std::vector<Photon> photons;

    double plasmaInternalEnergy;
    double photonPerMass; // := photon / proton * me / mp / (1 + Z_pm * me / mp)
    double meanParticleMass; // := mbar / me, where mbar := (np mp + ne me) / (np + ne) ~ mp
    double photonMeanFreePath; // with respect to eddy scale
    double zPlusMinus; // num (+ or -) over num protons
    double photonPerProton;
    double protonToElectronMassRatio;
};
