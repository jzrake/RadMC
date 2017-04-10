#ifndef RichardsonCascade_hpp
#define RichardsonCascade_hpp

#include "TabulatedFunction.hpp"



class RichardsonCascade
{
public:
    struct TimeScales
    {
        double getShortest();
        double eddyTime;
        double viscousTime;
        double comptonTime;
    };

    RichardsonCascade();
    ~RichardsonCascade();

    void advance (double dt);
    double getShortestTimeScale();
    TimeScales getTimeScales (int binIndex);

    TabulatedFunction spectralEnergy;

private:
    double cascadePower;
    double meanFreePath;
    double photonViscosity;
    double molecularViscosity;
};


#endif
