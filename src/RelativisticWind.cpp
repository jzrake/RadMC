#include <iostream>
#include <cmath>
#include "RelativisticWind.hpp"
// #define ADIABATIC




RelativisticWind::WindState::WindState (const RelativisticWind& wind, double r, double u, double s) :
r(r),
u(u),
specificWindPower (wind.specificWindPower),
initialFourVelocity (wind.initialFourVelocity),
adiabaticIndex (wind.adiabaticIndex)
{
    const double gm = wind.adiabaticIndex;
    const double e0 = wind.specificWindPower;
    const double u0 = wind.initialFourVelocity;
    const double r0 = 1.0;
    const double f0 = 1.0;
    const double g0 = std::sqrt (1 + u0 * u0);
    const double w0 = e0 / g0;
    const double m0 = w0 - 1;
    const double d0 = f0 / (r0 * r0 * u0);
    const double p0 = m0 * d0 * (gm - 1) / gm;
    const double s0 = s == 0.0 ? p0 / std::pow (d0, gm) : s;
    WindState::s = s0;

    d = f0 / (r * r * u);
    p = std::pow (d, gm) * s0;
    m = p / d * gm / (gm - 1);
    g = std::sqrt (1 + u * u);

    if (m < 0.0)
    {
        throw std::runtime_error (
            "Wind solution has negative enthalpy, mu = "
            + std::to_string (m)
            + " at r = "
            + std::to_string (r));
    }
}

double RelativisticWind::WindState::temperature() const
{
    const double gm = adiabaticIndex;
    const double Z = leptonsPerBaryon;
    const double X = photonsPerBaryon;
    const double T = (Z + P.mp / P.me) / (1 + Z + X) * m * (gm - 1) / gm;
    return T;
}

double RelativisticWind::WindState::properNumberDensity (WindState::Species species) const
{
    const double R = innerRadiusCm;
    const double F = luminosityPerSteradian / specificWindPower; // L = F eta (F is in erg / s)
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

double RelativisticWind::WindState::blackbodyPhotonsPerProton() const
{
    const auto P = PhysicsConstants();
    const double Z = leptonsPerBaryon;
    const double X = photonsPerBaryon;
    const double R = innerRadiusCm;
    const double F = luminosityPerSteradian / specificWindPower; // L = F eta (F is in erg / s)
    const double E0 = F / R / R / P.c; // base units of energy density (erg / cm^3)
    const double ug = 3 * p * E0 * X / (1 + Z + X); // radiative part of internal energy
    const double kT = std::pow (ug / (8. / 15 * std::pow (P.pi, 5.)) * std::pow (P.h * P.c, 3.), 1. / 4);
    const double zt = 1.202056903; // zeta(3)
    const double np = properNumberDensity (Species::baryon);
    const double ng = 16. * P.pi * zt * std::pow (kT / P.h / P.c, 3.);
    return ng / np;
}

double RelativisticWind::WindState::thomsonMeanFreePath (UnitVector nhat) const
{
    const double cosTheta = nhat.pitchAngleWith (propagationAngle);
    const double ne = properNumberDensity (Species::electron);
    const double dl = 1.0;
    const double dt = dl * (1 - u / g * cosTheta) * g * ne * P.st;
    return dl / dt;
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

void RelativisticWind::setEntropyProductionRate (double zeta)
{
    entropyProductionRate = zeta;
    resetSolverFunction();
}

RelativisticWind::WindState RelativisticWind::integrate (double outerRadius) const
{
#ifdef ADIABATIC
    solver.setValues (initialFourVelocity, 1.0);
    solver.integrate (outerRadius);

    const double r = solver.getT();
    const double u = solver.getY();

    return WindState (*this, r, u);
#else
    const double initialRadius = 1.0;
    const double initialEntropy = WindState (*this, initialRadius, initialFourVelocity).s;

    vsolver.setValues ({initialFourVelocity, initialEntropy}, initialRadius);
    vsolver.integrate (outerRadius);

    const double r = vsolver.getT();
    const double u = vsolver.getY()[0];
    const double s = vsolver.getY()[1];

    return WindState (*this, r, u, s);
#endif
}

std::vector<RelativisticWind::WindState> RelativisticWind::integrateTable (std::vector<double> radius) const
{
#ifdef ADIABATIC
    solver.setValues (initialFourVelocity, 1.0);
    solver.integrate (radius.empty() ? 0.0 : radius[0]);

    auto solution = std::vector<WindState>();

    for (int i = 0; i < radius.size(); ++i)
    {
        solver.integrate (radius[i]);

        const double r = solver.getT();
        const double u = solver.getY();
        const auto S = WindState (*this, r, u);
        solution.push_back(S);
    }
    return solution;
#else
    const double initialRadius = 1.0;
    const double initialEntropy = WindState (*this, initialRadius, initialFourVelocity).s;

    vsolver.setValues ({initialFourVelocity, initialEntropy}, initialRadius);
    vsolver.integrate (radius.empty() ? 0.0 : radius[0]);

    auto solution = std::vector<WindState>();

    for (int i = 0; i < radius.size(); ++i)
    {
        vsolver.integrate (radius[i]);

        const double r = vsolver.getT();
        const double u = vsolver.getY()[0];
        const double s = vsolver.getY()[1];
        const auto S = WindState (*this, r, u, s);
        solution.push_back(S);
    }
    return solution;
#endif
}

void RelativisticWind::resetSolverFunction()
{
#ifdef ADIABATIC
    const double gm = adiabaticIndex;
    const double e0 = specificWindPower;
    const double u0 = initialFourVelocity;
    const double r0 = 1.0;
    const double f0 = 1.0;

    const double g0 = std::sqrt (1 + u0 * u0);
    const double w0 = e0 / g0;
    const double m0 = w0 - 1;
    const double d0 = f0 / (r0 * r0 * u0);
    const double p0 = m0 * d0 * (gm - 1) / gm;
    const double s0 = p0 / std::pow (d0, gm);

    const double dsdlogr = 0.0;

    auto udot = [=] (double u, double r)
    {
        const double d = f0 / (r * r * u);
        const double p = std::pow (d, gm) * s0;
        const double g = std::sqrt (1 + u * u);
        const double m = p / d * gm / (gm - 1);
        const double w = 1 + m;

        const double v2 = u * u / g / g;
        const double a2 = gm * p / (d * w);

        const double P = (a2 - v2) * e0 / (gm - 1);
        const double Q = 2 * g * m;
        const double R = -(gm - 2) / (gm - 1) * g * m / gm * dsdlogr * u / r;
        const double du = -(R + Q * u / r) / P;

        return du;
    };
    solver.setTolerance (1e-10);
    solver.setFunction (udot);
#else
    const double gm = adiabaticIndex;
    const double e0 = specificWindPower;
    const double f0 = 1.0;
    const double dsdlogr = entropyProductionRate;

    auto udot = [=] (std::valarray<double> y, double r) -> std::valarray<double>
    {
        const double u = y[0];
        const double s = y[1];
        const double d = f0 / (r * r * u);
        const double p = std::pow (d, gm) * s;
        const double g = std::sqrt (1 + u * u);
        const double m = p / d * gm / (gm - 1);
        const double w = 1 + m;

        const double v2 = u * u / g / g;
        const double a2 = gm * p / (d * w);

        const double P = (a2 - v2) * e0 / (gm - 1);
        const double Q = 2 * g * m;
        const double R = -(gm - 2) / (gm - 1) * g * m / gm * dsdlogr * u / r;
        const double du = -(R + Q * u / r) / P;
        const double ds = dsdlogr / r;

        return {du, ds};
    };

    vsolver.setTolerance (1e-10);
    vsolver.setFunction (udot);
#endif
}
