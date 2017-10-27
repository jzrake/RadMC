#include <cmath>
#include "RelativisticWind.hpp"




RelativisticWind::WindState::WindState (const RelativisticWind& wind, double r, double u) : r(r), u(u)
{
    const double gm = wind.adiabaticIndex;
    const double e0 = wind.windPower;
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

RelativisticWind::RelativisticWind()
{
    resetSolverFunction();
}

void RelativisticWind::setWindPower (double eta)
{
    windPower = eta;
    resetSolverFunction();
}

void RelativisticWind::setInitialFourVelocity (double u0)
{
    initialFourVelocity = u0;
    resetSolverFunction();
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
    const double e0 = windPower;
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
