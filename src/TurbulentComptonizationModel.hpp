#pragma once

#include "TabulatedFunction.hpp"
#include "FourVector.hpp"
#include "RichardsonCascade.hpp"




class TurbulentComptonizationModel
{
public:




    struct Config
    {
        std::string outdir = ".";
        double tmax = 1.0;
        double cpi = 1.0; // Checkpoint interval
        double tsi = 1.0; // Time series interval
        double theta = 0.01; // electron temperature in units of electron rest mass
        double ell_star = 1e-3; // photon mean free path
        double nphot = 4; // log10 of photon number
        double nphot_per_mass = 1e-1; // Photon number per unit mass (1 means mp / me photons per baryon)
        double nelec = 5; // log10 of electron number (currently for diagnostics only)
        double ephot = 1e-3;
        double beta_turb = 0.5;
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

    double getTimestep() const;
    void advance (double dt);


    std::vector<double> getCascadeWaveNumbers() const;
    std::vector<double> getCascadePowerSpectrum() const;
    std::vector<Photon> getPhotons() const;
    double getFluidVelocityAtPhotonMeanFreePathScale() const;
    double getTotalPhotonEnergyMC() const;
    double getMeanPhotonEnergy() const;    
    double getSpecificPhotonEnergy() const;
    double getElectronTemperature() const;
    double getPhotonTemperature() const;
    double getSpecificInternalEnergy (double electronTemperature) const;
    double getComptonCoolingTime() const;

private:
    Electron sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta) const;
    Electron sampleElectronForScatteringInParcel (const Photon& photon, RandomVariable& electronGammaBeta) const;
    FourVector doComptonScattering (Photon& photon, Electron& electron) const;
    void computeNextPhotonScatteringAndParcelVelocity (Photon& photon) const;
    void regenerateElectronVelocityDistribution (double electronTemperature);

    Config config;

    RandomVariable electronGammaBeta;
    RandomVariable photonEnergy;
    std::vector<Photon> photons;

    double plasmaInternalEnergy;
    double fluidKineticEnergy;
    double photonPerMass; // := photon / proton * me / mp / (1 + Z_pm * me / mp)
    double meanParticleMass; // := mbar / me, where mbar := (np mp + ne me) / (np + ne) ~ mp
    double photonMeanFreePath; // with respect to eddy scale
    double zPlusMinus; // num (+ or -) over num protons
    double photonPerProton;
    double protonToElectronMassRatio;
    
    RichardsonCascade cascadeModel;
    bool coldElectrons;
    double simulationTime;
};

