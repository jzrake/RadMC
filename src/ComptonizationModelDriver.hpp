#ifndef ComptonizationModelDriver_hpp
#define ComptonizationModelDriver_hpp

#include "SimulationDriver.hpp"
#include "TabulatedFunction.hpp"
#include "FourVector.hpp"



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




class ComptonizationModelDriver : public SimulationDriver
{
public:
    ComptonizationModelDriver();
    void makeUserParameters (Variant::NamedValues& params) override;
    void configureFromParameters() override;
    double getTimestep() const override;
    bool shouldContinue() const override;
    void advance (double dt) override;
    bool shouldRecordIterationInTimeSeries() const override;
    double getRecordForTimeSeries (std::string) const override;
    bool shouldWriteOutput() const override;
    void writeOutput() const override;

private:
    Electron sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta) const;
    Electron sampleElectronForScatteringInParcel (const Photon& photon, RandomVariable& electronGammaBeta) const;
    void computeNextPhotonScatteringAndParcelVelocity (Photon& photon) const;
    FourVector doComptonScattering (Photon& photon, Electron& electron) const;
    double getMeanPhotonEnergy() const;

    RandomVariable electronGammaBeta;
    RandomVariable photonEnergy;
    std::vector<Photon> photons;
    Electron electronPopulation;
    double fluidKineticEnergy;
};

#endif