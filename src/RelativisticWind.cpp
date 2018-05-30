#include <iostream>
#include <cmath>
#include "RelativisticWind.hpp"




// ============================================================================
RelativisticWind::WindState::WindState (double e, double u, double r, double heatingRate)
: e(e)
, u(u)
, r(r)
, heatingRate (heatingRate)
{
    const double gm = 4. / 3;
    const double f0 = 1.0;  

    g = std::sqrt (1 + u * u);
    h = e / g;
    w = h - 1;
    d = f0 / (r * r * u);
    p = d * w * (gm - 1) / gm;
    s = std::log (p / std::pow (d, gm));

    const double v2 = u * u / g / g;
    const double a2 = gm * p / (d * h);
    M2 = v2 / a2;

    if (w < 0.0)
    {
        throw std::runtime_error (
            "Wind solution has negative enthalpy, mu = "
            + std::to_string(w)
            + " at r = "
            + std::to_string(r));
    }
}

double RelativisticWind::WindState::radius() const
{
    return r * innerRadiusCm;
}

double RelativisticWind::WindState::dLogedLogr() const
{
    return (w / h) * heatingRateXi();
}

double RelativisticWind::WindState::dLogudLogr() const
{
    return (2 - heatingRateXi()) / (M2 - 1);
}

double RelativisticWind::WindState::dLogwdLogr() const
{
    return -M2 / 3 * dLogudLogr() + heatingRateXi();
}

double RelativisticWind::WindState::heatingRateXi() const
{
    // For this formula, heatingRate is Chris's tilde-xi, and must be multiplied by eta0.
    return heatingRate / (u * w);
}

double RelativisticWind::WindState::properNumberDensity (Species species) const
{
    const double R = innerRadiusCm;
    const double F = luminosityPerSteradian / e; // L = F eta (F is in erg / s)
    const double Z = leptonsPerBaryon;
    const double D = d * F / R / R / P.c; // rest-mass density (erg / cm^3)
    const double np = D / (P.mp * P.c * P.c * (1 + Z * P.me / P.mp));

    switch (species)
    {
        case Species::baryon: return np;
        case Species::electron: return np * leptonsPerBaryon;
        case Species::photon: return np * photonsPerBaryon;
    }
}

double RelativisticWind::WindState::jetOpticalDepth() const
{
    const double ne = properNumberDensity (Species::electron);
    return ne * P.st * radius() / g;
}

double RelativisticWind::WindState::blackbodyPhotonsPerProton() const
{
    const double Z = leptonsPerBaryon;
    const double X = photonsPerBaryon; // CHECK: why does this function use photonsPerBaryon?
    const double R = innerRadiusCm;
    const double F = luminosityPerSteradian / e; // L = F eta (F is in erg / s)
    const double E0 = F / R / R / P.c; // base units of energy density (erg / cm^3)
    const double ug = 3 * p * E0 * X / (1 + Z + X); // radiative part of internal energy
    const double kT = std::pow (ug / (8. / 15 * std::pow (P.pi, 5.)) * std::pow (P.h * P.c, 3.), 1. / 4);
    const double zt = 1.202056903; // zeta(3)
    const double np = properNumberDensity (Species::baryon);
    const double ng = 16. * P.pi * zt * std::pow (kT / P.h / P.c, 3.);
    return ng / np;
}

double RelativisticWind::WindState::photonTemperature() const
{
    const double gm = 4. / 3;
    const double Z = leptonsPerBaryon;
    const double X = photonsPerBaryon;
    const double T = (Z + P.mp / P.me) / (1 + Z + X) * w * (gm - 1) / gm;
    return T;
}

double RelativisticWind::WindState::comptonParameter() const
{
    const double alpha = 2 + dLogudLogr();
    return dLogwdLogr() + alpha / 3.;
}

double RelativisticWind::WindState::deltaTheta() const
{
    return comptonParameter() / jetOpticalDepth() / 4.;
}

double RelativisticWind::WindState::electronTemperature() const
{
    return photonTemperature() + deltaTheta();
}

double RelativisticWind::WindState::thomsonMeanFreePath (UnitVector nhat) const
{
    const double cosTheta = nhat.pitchAngleWith (propagationAngle);
    const double ne = properNumberDensity (Species::electron);
    const double dl = 1.0;
    const double dt = dl * (1 - u / g * cosTheta) * g * ne * P.st;
    return dl / dt;
}

double RelativisticWind::WindState::thomsonMeanFreePathComoving() const
{
    const double ne = properNumberDensity (Species::electron);
    return 1.0 / (ne * P.st);
}

double RelativisticWind::WindState::radiationViscosity() const
{
    return w / (1 + w) * thomsonMeanFreePathComoving() * P.c;
}

double RelativisticWind::WindState::causallyConnectedScale() const
{
    return innerRadiusCm * r / g;
}

FourVector RelativisticWind::WindState::fourVelocity() const
{
    return FourVector::fromGammaBetaAndUnitVector (u, propagationAngle);
}




// ============================================================================
RelativisticWind::RelativisticWind()
{
    resetSolverFunction();
}

void RelativisticWind::setSpecificWindPower (double eta)
{
    specificWindPower = eta;
    resetSolverFunction();
}

void RelativisticWind::setInitialFourVelocity (double u0)
{
    initialFourVelocity = u0;
    resetSolverFunction();
}

void RelativisticWind::setHeatingRate (double zeta)
{
    heatingRate = zeta;
    resetSolverFunction();
}

RelativisticWind::WindState RelativisticWind::integrate (double outerRadius) const
{
    const double initialRadius = 1.0;

    vsolver.setValues ({specificWindPower, initialFourVelocity}, initialRadius);
    vsolver.integrate (outerRadius);
    const double e = vsolver.getY()[0]; // eta
    const double u = vsolver.getY()[1]; // four-velocity
    const double r = vsolver.getT();    // radius
    return WindState (e, u, r, heatingRate);
}

std::vector<RelativisticWind::WindState> RelativisticWind::integrateRange (double outerRadius) const
{
    const double r0 = 1.0;
    const double r1 = outerRadius;
    const int numberOfBins = 400;
    std::vector<double> r;

    for (int i = 0; i < numberOfBins; ++i)
    {
        r.push_back (r0 * std::pow (r1 / r0, double(i) / (numberOfBins - 1)));
    }
    return integrateTable(r);
}

std::vector<RelativisticWind::WindState> RelativisticWind::integrateTable (std::vector<double> radius) const
{
    const double initialRadius = 1.0;

    vsolver.setValues ({specificWindPower, initialFourVelocity}, initialRadius);
    vsolver.integrate (radius.empty() ? 0.0 : radius[0]);

    auto solution = std::vector<WindState>();

    for (int i = 0; i < radius.size(); ++i)
    {
        vsolver.integrate (radius[i]);
        const double e = vsolver.getY()[0]; // eta
        const double u = vsolver.getY()[1]; // four-velocity
        const double r = vsolver.getT();    // radius
        const auto S = WindState (e, u, r, heatingRate);
        solution.push_back(S);
    }
    return solution;
}

void RelativisticWind::resetSolverFunction()
{
    auto udot = [=] (std::valarray<double> y, double r) -> std::valarray<double>
    {
        const double e = y[0];
        const double u = y[1];
        WindState S (e, u, r, heatingRate);
        return {S.dLogedLogr() * S.e / S.r, S.dLogudLogr() * S.u / S.r};
    };
    vsolver.setTolerance (1e-8);
    vsolver.setFunction (udot);
}
