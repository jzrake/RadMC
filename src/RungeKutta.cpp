#include <cmath>
#include "RungeKutta.hpp"




RungeKutta::RungeKutta (DerivativeFunction ydot) : ydot (ydot)
{

}

void RungeKutta::setFunction (DerivativeFunction ydotToUse)
{
    ydot = ydotToUse;
}

void RungeKutta::setValues (double y, double t)
{
    yn = y;
    tn = t;
}

void RungeKutta::setTolerance (double relativeErrorToUse)
{
    relativeError = relativeErrorToUse;
}

double RungeKutta::getT() const
{
    return tn;
}

double RungeKutta::getY() const
{
    return yn;
}

double RungeKutta::integrate (double t)
{
    if (ydot == nullptr)
    {
        throw std::logic_error ("RungeKutta: no ydot function has been set");
    }

    divisions = 1;

    double dt = t - tn;
    double y1 = integrateSubdivided (yn, tn, dt, divisions);

    while (divisions < (1 << 20))
    {
        double y2 = integrateSubdivided (yn, tn, dt, divisions * 2);
        double error = std::fabs (y2 - y1) / std::fabs (y2);

        if (error < relativeError)
        {
            tn = t;
            yn = y2;
            return yn;        
        }
        y1 = y2;
        divisions *= 2;
    }
    throw std::runtime_error ("RungeKutta::integrate is taking too long to converge");
}

int RungeKutta::getDivisionsForLast() const
{
    return divisions;
}

double RungeKutta::evaluateDeltaY (double y0, double t0, double dt) const
{
    const double k0 = 0.0;
    const double k1 = ydot (y0 + k0 * 0.0 * dt, t0 + 0.0 * dt);
    const double k2 = ydot (y0 + k1 * 0.5 * dt, t0 + 0.5 * dt);
    const double k3 = ydot (y0 + k2 * 0.5 * dt, t0 + 0.5 * dt);
    const double k4 = ydot (y0 + k3 * 1.0 * dt, t0 + 1.0 * dt);

    return (k1 + k2 * 2 + k3 * 2 + k4) * dt / 6;
}

double RungeKutta::integrateSubdivided (double y0, double t0, double dt, int steps) const
{
    for (int n = 0; n < steps; ++n)
    {
        y0 += evaluateDeltaY (y0, t0, dt / steps);
        t0 += dt / steps;
    }
    return y0;
}
