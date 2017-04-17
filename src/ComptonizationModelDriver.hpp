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

    /**
    Return the electron four-velocity (this is the same as its momentum since
    electron mass is 1).
    */
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
    void doComptonScattering (Photon& photon, Electron& electron);
    double getMeanPhotonEnergy() const;

    RandomVariable electronGammaBeta;
    RandomVariable photonEnergy;
    std::vector<Photon> photons;
};

#endif