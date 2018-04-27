#include <iostream>
#include <cmath>
#include "RelativisticWind.hpp"
// #define ADIABATIC




RelativisticWind::WindState::WindState (const RelativisticWind& wind, double r, double u, double ef)
: r(r)
, u(u)
, specificWindPower (wind.specificWindPower)
, initialFourVelocity (wind.initialFourVelocity)
, adiabaticIndex (wind.adiabaticIndex)
{
    const double gm = wind.adiabaticIndex;
    const double e0 = wind.specificWindPower;
    const double f0 = 1.0;  

    g = std::sqrt (1 + u * u);

    h = e0 / g;
    n = ef / g;
    w = h - n - 1;
    p = d * w * (gm - 1) / gm;
    e = ef;
    d = f0 / (r * r * u);
    p = d * w * (gm - 1) / gm;
    s = std::log (p / std::pow (d, gm));

    if (w < 0.0)
    {
        throw std::runtime_error (
            "Wind solution has negative enthalpy, mu = "
            + std::to_string(w)
            + " at r = "
            + std::to_string(r));
    }
}

double RelativisticWind::WindState::temperature() const
{
    const double gm = adiabaticIndex;
    const double Z = leptonsPerBaryon;
    const double X = photonsPerBaryon;
    const double T = (Z + P.mp / P.me) / (1 + Z + X) * w * (gm - 1) / gm;
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
    const double X = photonsPerBaryon; // CHECK: why does this function use photonsPerBaryon?
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

double RelativisticWind::WindState::thomsonMeanFreePathComoving() const
{
    const double ne = properNumberDensity (Species::electron);
    return 1.0 / (ne * P.st);
}

double RelativisticWind::WindState::radiationViscosity() const
{
    const auto P = PhysicsConstants();
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

void RelativisticWind::setSpecificWindPower (double eta, double etaFree)
{
    specificWindPower = eta;
    specificFreePower = etaFree;
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
#ifdef ADIABATIC
    solver.setValues (initialFourVelocity, 1.0);
    solver.integrate (outerRadius);

    const double r = solver.getT();
    const double u = solver.getY();

    return WindState (*this, r, u);
#else
    const double initialRadius = 1.0;

    vsolver.setValues ({initialFourVelocity, specificFreePower}, initialRadius);
    vsolver.integrate (outerRadius);

    const double r = vsolver.getT();
    const double u = vsolver.getY()[0];
    const double f = vsolver.getY()[1];

    return WindState (*this, r, u, f);
#endif
}

std::vector<RelativisticWind::WindState> RelativisticWind::integrateRange (double outerRadius) const
{
    const double r0 = 1.0;
    const double r1 = outerRadius;
    const int numberOfBins = 100;
    std::vector<double> r;

    for (int i = 0; i < numberOfBins; ++i)
    {
        r.push_back (r0 * std::pow (r1 / r0, double(i) / (numberOfBins - 1)));
    }
    return integrateTable(r);
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

    vsolver.setValues ({initialFourVelocity, specificFreePower}, initialRadius);
    vsolver.integrate (radius.empty() ? 0.0 : radius[0]);

    auto solution = std::vector<WindState>();

    for (int i = 0; i < radius.size(); ++i)
    {
        vsolver.integrate (radius[i]);

        const double r = vsolver.getT();
        const double u = vsolver.getY()[0];
        const double e = vsolver.getY()[1]; // eta-free
        const auto S = WindState (*this, r, u, e);
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
    const double f0 = 1.0;
    const double e0 = specificWindPower;
    const double gm = adiabaticIndex;
    const double dlogedlogr = -heatingRate;

    auto udot = [=] (std::valarray<double> y, double r) -> std::valarray<double>
    {
        const double u = y[0];
        const double e = y[1];
        const double d = f0 / (r * r * u);
        const double g = std::sqrt (1 + u * u);
        const double h = e0 / g;
        const double n = e / g;
        const double w = h - n - 1;
        const double p = d * w * (gm - 1) / gm;

        const double v2 = u * u / g / g;
        const double a2 = gm * p / (d * h);
        const double M2 = v2 / a2;

        const double dlogudlogr = (2 + n / w * dlogedlogr) / (M2 - 1) / (1 + v2 * n / w / (M2 - 1));
        const double dudr = dlogudlogr * u / r;
        const double dedr = dlogedlogr * e / r;

        // std::cout << dlogedlogr << " " << e << std::endl;

        return {dudr, dedr};
    };

    vsolver.setTolerance (1e-8);
    vsolver.setFunction (udot);
#endif
}
