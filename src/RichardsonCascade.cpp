//#include <iostream>
#include <cmath>
#include "RichardsonCascade.hpp"




// ============================================================================
double RichardsonCascade::TimeScales::getShortest()
{
    double timeScales[] = {eddyTurnoverTime, viscousDampingTime, comptonDragTime};
    return *std::min_element (timeScales, timeScales + 3);
}




// ============================================================================
RichardsonCascade::RichardsonCascade() :
spectralEnergy (1, 100000, 128, TabulatedFunction::useEqualBinWidthsLogarithmic)
{
    spectralEnergy[0] = 1;
    cascadePower = 1;
    photonMeanFreePath = 1e-3;
    radiativeEnergyDensity = 8e-5;
}

RichardsonCascade::~RichardsonCascade()
{

}

void RichardsonCascade::advance (double dt)
{
    std::vector<double> energyFlux;
    std::vector<double> dampingLoss;

    spectralEnergy[0] += cascadePower * dt;

    for (int n = 0; n < spectralEnergy.size(); ++n)
    {
        double dE = spectralEnergy[n];
        TimeScales T = getTimeScales (n);

        energyFlux.push_back (dE / T.eddyTurnoverTime);
        dampingLoss.push_back (dE * (1 / T.viscousDampingTime + 1 / T.comptonDragTime));
    }

    for (int n = 0; n < spectralEnergy.size() - 1; ++n)
    {
        // Conservative
        spectralEnergy[n]     -= dt * energyFlux[n];
        spectralEnergy[n + 1] += dt * energyFlux[n];

        // Non-conservative
        spectralEnergy[n] -= dt * dampingLoss[n];
    }
}

double RichardsonCascade::getShortestTimeScale() const
{
    double shortestTime = 0;
    bool first = true;

    for (int n = 0; n < spectralEnergy.size(); ++n)
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

RichardsonCascade::TimeScales RichardsonCascade::getTimeScales (int binIndex) const
{
    const double dE = spectralEnergy[binIndex];
    const double eddySpeed = std::sqrt (dE);
    const double eddyScale = 1.0 / spectralEnergy.getBinEdge (binIndex);
    const double viscosity = photonMeanFreePath * radiativeEnergyDensity; // Note: needs 8 / 27 in front

    TimeScales T;
    T.eddyTurnoverTime = eddyScale / eddySpeed;
    T.viscousDampingTime = eddyScale * eddyScale / viscosity;
    T.comptonDragTime = 3. / 4 * (photonMeanFreePath / radiativeEnergyDensity); // Note: needs a Zpm

    return T;
}

std::vector<double> RichardsonCascade::getEddyTurnoverTime() const
{
    std::vector<double> eddyTurnoverTime;

    for (int n = 0; n < spectralEnergy.size(); ++n)
    {
        TimeScales T = getTimeScales (n);
        eddyTurnoverTime.push_back (T.eddyTurnoverTime);
    }
    return eddyTurnoverTime;
}

std::vector<double> RichardsonCascade::getDampingTime() const
{
    std::vector<double> dampingTime;

    for (int n = 0; n < spectralEnergy.size(); ++n)
    {
        TimeScales T = getTimeScales (n);
        dampingTime.push_back (std::min (T.viscousDampingTime, T.comptonDragTime));
    }
    return dampingTime;
}

