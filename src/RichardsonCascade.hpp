#ifndef RichardsonCascade_hpp
#define RichardsonCascade_hpp

#include "TabulatedFunction.hpp"



class RichardsonCascade
{
public:
    struct TimeScales
    {
        double getShortest();
        double eddyTurnoverTime;
        double viscousDampingTime;
        double comptonDragTime;
    };

    RichardsonCascade();
    ~RichardsonCascade();

    void advance (double dt);
    double getShortestTimeScale() const;
    TimeScales getTimeScales (int binIndex) const;
    std::vector<double> getEddyTurnoverTime() const;
    std::vector<double> getDampingTime() const;

    TabulatedFunction spectralEnergy;

private:
    double cascadePower;           // should generally be 1
    double photonMeanFreePath;     // in units of the outer scale
    double radiativeEnergyDensity; // in units of rho c^2
};


#endif
