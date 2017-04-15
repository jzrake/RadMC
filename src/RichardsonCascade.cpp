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
powerSpectrum (1, 1e4, 256, TabulatedFunction::useEqualBinWidthsLogarithmic)
{
    cascadePower = 1;
    photonMeanFreePath = 1e-3;
    radiativeEnergyDensity = 1e-5;
}

RichardsonCascade::~RichardsonCascade()
{

}

void RichardsonCascade::advance (double dt)
{
    std::vector<double> energyFlux;

    energyFlux.push_back (cascadePower);

    for (int n = 0; n < powerSpectrum.size(); ++n)
    {
        const double kn = powerSpectrum.getBinEdge (n);
        const double Pk = powerSpectrum[n];
        const double Fk = std::pow (Pk, 3. / 2) * std::pow (kn, 5. / 2);
        energyFlux.push_back (Fk);
    }

    for (int n = 0; n < powerSpectrum.size() - 1; ++n)
    {
        const double kn = powerSpectrum.getBinEdge (n);
        const double dk = powerSpectrum.getBinWidth (n);
        const double Pk = powerSpectrum[n];
        const double nu = photonMeanFreePath * radiativeEnergyDensity;
        const double viscousPower = 2 * kn * kn * nu * Pk;

        powerSpectrum[n] -= dt * (energyFlux[n + 1] - energyFlux[n]) / dk;
        powerSpectrum[n] -= dt * viscousPower;
    }
}

RichardsonCascade::TimeScales RichardsonCascade::getTimeScales (int binIndex) const
{
    const double Pk = powerSpectrum[binIndex];
    const double el = 1.0 / powerSpectrum.getBinEdge (binIndex);
    const double vk = std::sqrt (Pk / el);
    const double nu = photonMeanFreePath * radiativeEnergyDensity; // Note: needs an 8 / 27

    TimeScales T;
    T.eddyTurnoverTime = el / vk;
    T.viscousDampingTime = el * el / nu;
    T.comptonDragTime = photonMeanFreePath / radiativeEnergyDensity; // Note: needs a Zpm and a 3/4
    T.effectiveDampingTime = el < photonMeanFreePath ? T.comptonDragTime : T.viscousDampingTime;

    return T;
}

double RichardsonCascade::getShortestTimeScale() const
{
    double shortestTime = 0;
    bool first = true;

    for (int n = 0; n < powerSpectrum.size(); ++n)
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

    for (int n = 0; n < powerSpectrum.size(); ++n)
    {
        TimeScales T = getTimeScales (n);
        eddyTurnoverTime.push_back (T.eddyTurnoverTime);
    }
    return eddyTurnoverTime;
}

std::vector<double> RichardsonCascade::getDampingTime() const
{
    std::vector<double> dampingTime;

    for (int n = 0; n < powerSpectrum.size(); ++n)
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
    double viscosity = photonMeanFreePath * radiativeEnergyDensity; // Note: needs an 8 / 27
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

double RichardsonCascade::getTotalEnergy() const
{
    double E = 0;

    for (int n = 0; n < powerSpectrum.size() - 1; ++n)
    {
        E += powerSpectrum[n] * powerSpectrum.getBinWidth (n);
    }
    return E;
}
