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
    Photon (FourVector momentum) : momentum (momentum) {}

    /**
    Update the position four-vector based on the momentum and time step dt.
    */
    void advancePosition (double dt)
    {
        UnitVector nhat = momentum.getUnitThreeVector();
        position += FourVector (dt, nhat.getX() * dt, nhat.getY() * dt, nhat.getZ() * dt);
    }

    FourVector position;
    FourVector momentum;
};




class ComptonizationModelDriver : public SimulationDriver
{
public:
    ComptonizationModelDriver();
    void makeUserParameters (Variant::NamedValues& params) override;
    void configureFromParameters() override;
    void printStartupMessage() const override;
    double getTimestep() const override;
    bool shouldContinue() const override;
    void advance (double dt) override;
    bool shouldRecordIterationInTimeSeries() const override;
    double getRecordForTimeSeries (std::string) const override;
    bool shouldWriteOutput() const override;
    void writeOutput () const override;

private:
    Electron sampleElectronForScattering (const Photon& photon, RandomVariable& electronGammaBeta);
    Electron sampleElectronForScattering (const Photon& photon, const FourVector::Field& field);
    void doComptonScattering (Photon& photon, Electron& electron);
    double getMeanPhotonEnergy() const;

    RandomVariable electronGammaBeta;
    RandomVariable photonEnergy;
    std::vector<Photon> photons;
};

#endif