#include <cmath>
#include "RelativisticWind.hpp"




RelativisticWind::WindState::WindState (const RelativisticWind& wind, double r, double u) :
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
    const double s0 = p0 / std::pow (d0, gm);

    d = f0 / (r * r * u);
    m = p / d * gm / (gm - 1);
    p = std::pow (d, gm) * s0;
    g = std::sqrt (1 + u * u);
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
    const double F = luminosityPerSteradian / specificWindPower; // L = 4 pi F eta (F is in erg / s)
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

RelativisticWind& RelativisticWind::setSpecificWindPower (double eta)
{
    specificWindPower = eta;
    resetSolverFunction();
    return *this;
}

RelativisticWind& RelativisticWind::setInitialFourVelocity (double u0)
{
    initialFourVelocity = u0;
    resetSolverFunction();
    return *this;
}

RelativisticWind::WindState RelativisticWind::integrate (double outerRadius) const
{
    solver.setValues (initialFourVelocity, 1.0);
    solver.integrate (outerRadius);

    const double r = solver.getT();
    const double u = solver.getY();

    return WindState (*this, r, u);
}

std::vector<RelativisticWind::WindState> RelativisticWind::integrate (std::vector<double> radius) const
{
    if (radius.empty())
    {
        return std::vector<WindState>();
    }
    auto solution = std::vector<WindState>();

    solver.setValues (initialFourVelocity, 1.0);
    solver.integrate (radius[0]);

    for (int i = 0; i < radius.size(); ++i)
    {
        solver.integrate (radius[i]);
        const double r = solver.getT();
        const double u = solver.getY();
        solution.push_back (WindState (*this, r, u));
    }
    return solution;
}

void RelativisticWind::resetSolverFunction()
{
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
        const double R = 0.0;
        const double du = -(R + Q * u / r) / P;

        return du;
    };
    solver.setFunction (udot);
}
