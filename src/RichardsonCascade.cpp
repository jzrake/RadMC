//#include <iostream>
#include <cmath>
#include "RichardsonCascade.hpp"




// ============================================================================
double RichardsonCascade::TimeScales::getShortest()
{
    return viscousTime < eddyTime ? viscousTime : eddyTime;
}




// ============================================================================
RichardsonCascade::RichardsonCascade() :
spectralEnergy (1, 10000, 128, TabulatedFunction::useEqualBinWidthsLogarithmic)
{
    spectralEnergy[0] = 1;
    cascadePower = 1;
    meanFreePath = 1e-3;
    photonViscosity = 2e-5;
    molecularViscosity = 1e-9;
}

RichardsonCascade::~RichardsonCascade()
{

}

void RichardsonCascade::advance (double dt)
{
    std::vector<double> energyFlux;
    std::vector<double> viscousLoss;

    spectralEnergy[0] += cascadePower * dt;

    for (int n = 0; n < spectralEnergy.size() - 1; ++n)
    {
        double dE = spectralEnergy[n];
        TimeScales T = getTimeScales (n);

        energyFlux.push_back (dE / T.eddyTime);
        viscousLoss.push_back (dE / T.viscousTime);
    }

    for (int n = 0; n < energyFlux.size(); ++n)
    {
        // Conservative
        spectralEnergy[n]     -= dt * energyFlux[n];
        spectralEnergy[n + 1] += dt * energyFlux[n];

        // Non-conservative
        spectralEnergy[n] -= dt * viscousLoss[n];
    }
}

double RichardsonCascade::getShortestTimeScale()
{
    double shortestTime = 0;
    bool first = true;

    for (int n = 0; n < spectralEnergy.size() - 1; ++n)
    {
        TimeScales T = getTimeScales (n);

        if (first || T.getShortest() < shortestTime)
        {
            first = false;
            shortestTime = T.getShortest();
        }
    }
    return shortestTime;
}

RichardsonCascade::TimeScales RichardsonCascade::getTimeScales (int binIndex)
{
    double dE = spectralEnergy[binIndex];
    double eddySpeed = std::sqrt (dE);
    double eddyScale = 1.0 / spectralEnergy.getBinEdge (binIndex);
    double eddyTime = eddyScale / eddySpeed;
    double viscosity = eddyScale < meanFreePath ? molecularViscosity : photonViscosity;
    double viscousTime = eddyScale * eddyScale / viscosity;

    TimeScales T;
    T.eddyTime = eddyTime;
    T.viscousTime = viscousTime;
    T.comptonTime = 0;

    return T;
}
