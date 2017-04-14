#include <iostream>
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
spectralEnergy (1, 1e7, 256, TabulatedFunction::useEqualBinWidthsLogarithmic)
{
    spectralEnergy[0] = 1;
    cascadePower = 1;
    photonMeanFreePath = 1e-3;
    radiativeEnergyDensity = 1e-5;
}

RichardsonCascade::~RichardsonCascade()
{

}

void RichardsonCascade::advance (double dt)
{
    // There are N energy bins
    std::vector<double> energyFlux; // N + 1 fluxes (the zero flux is the power)
    std::vector<double> dampingLoss; // N damping losses

    energyFlux.push_back (cascadePower);

    for (int n = 0; n < spectralEnergy.size(); ++n)
    {
        double dE = spectralEnergy[n];
        TimeScales T = getTimeScales (n);

        energyFlux.push_back (dE / T.eddyTurnoverTime);
        dampingLoss.push_back (dE / T.effectiveDampingTime);
    }

    for (int n = 0; n < spectralEnergy.size(); ++n)
    {
        // Conservative
        spectralEnergy[n] -= dt * (energyFlux[n + 1] - energyFlux[n]);

        // Non-conservative
        spectralEnergy[n] -= dt * dampingLoss[n];
    }
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
    T.comptonDragTime = photonMeanFreePath / radiativeEnergyDensity; // Note: needs a Zpm and a 3/4
    T.effectiveDampingTime = eddyScale < photonMeanFreePath ? T.comptonDragTime : T.viscousDampingTime;

    return T;
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

double RichardsonCascade::getPhotonMeanFreePathScale() const
{
    return photonMeanFreePath;
}

double RichardsonCascade::getFiducialViscousScale() const
{
    double outerScale = 1.0;
    double rmsEddySpeed = 1.0;
    double viscosity = photonMeanFreePath * radiativeEnergyDensity; // Note: needs 8 / 27 in front
    double Re = outerScale * rmsEddySpeed / viscosity;
    return outerScale * std::pow (Re, -3. / 4);
}

double RichardsonCascade::getFiducialComptonScale() const
{
    double outerScale = 1.0;
    double rmsEddySpeed = 1.0;
    double t0 = outerScale / rmsEddySpeed;
    double tc = photonMeanFreePath / radiativeEnergyDensity; // Note: needs a Zpm and a 3/4
    return outerScale * std::pow (tc / t0, 3. / 2);
}

double RichardsonCascade::getFiducialComptonPower() const
{
    double outerScale = 1.0;
    double rmsEddySpeed = 1.0;
    double tc = photonMeanFreePath / radiativeEnergyDensity; // Note: needs a Zpm and a 3/4
    double vs = rmsEddySpeed * std::pow (photonMeanFreePath / outerScale, 1. / 3);
    return (vs * vs / tc) / cascadePower;
}
