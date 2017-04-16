#ifndef TurbulenceModelDriver_hpp
#define TurbulenceModelDriver_hpp

#include "SimulationDriver.hpp"
#include "RichardsonCascade.hpp"




class TurbulenceModelDriver : public SimulationDriver
{
public:
    void makeUserParameters (Variant::NamedValues& params) override;
    void configureFromParameters() override;
    void printStartupMessage() const override;
    double getTimestep() const override;
    void advance (double dt) override;
    bool shouldContinue() const override;
    bool shouldWriteOutput() const override;
    void writeOutput (std::string filename) const override;
    std::string makeOutputFilename() const override;

private:
    RichardsonCascade cascade;
};

#endif