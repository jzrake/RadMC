#include <cmath>
#include "RungeKutta.hpp"




// ============================================================================
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




// ============================================================================
RungeKuttaVector::RungeKuttaVector (DerivativeFunction ydot) : ydot (ydot)
{

}

void RungeKuttaVector::setFunction (DerivativeFunction ydotToUse)
{
    ydot = ydotToUse;
}

void RungeKuttaVector::setValues (std::valarray<double> y, double t)
{
    yn = y;
    tn = t;
}

void RungeKuttaVector::setTolerance (double relativeErrorToUse)
{
    relativeError = relativeErrorToUse;
}

double RungeKuttaVector::getT() const
{
    return tn;
}

std::valarray<double> RungeKuttaVector::getY() const
{
    return yn;
}

std::valarray<double> RungeKuttaVector::integrate (double t)
{
    if (ydot == nullptr)
    {
        throw std::logic_error ("RungeKuttaVector: no ydot function has been set");
    }

    divisions = 1;

    auto dt = t - tn;
    auto y1 = integrateSubdivided (yn, tn, dt, divisions);

    while (divisions < (1 << 20))
    {
        auto y2 = integrateSubdivided (yn, tn, dt, divisions * 2);
        auto error = std::abs (y2 - y1); // / std::abs (y2);

        if (error.sum() < relativeError)
        {
            tn = t;
            yn = y2;
            return yn;        
        }
        y1 = y2;
        divisions *= 2;
    }
    throw std::runtime_error ("RungeKuttaVector::integrate is taking too long to converge");
}

int RungeKuttaVector::getDivisionsForLast() const
{
    return divisions;
}

std::valarray<double> RungeKuttaVector::evaluateDeltaY (std::valarray<double> y0, double t0, double dt) const
{
    const auto k0 = std::valarray<double>(0.0, y0.size());
    const auto k1 = ydot (y0 + k0 * 0.0 * dt, t0 + 0.0 * dt);
    const auto k2 = ydot (y0 + k1 * 0.5 * dt, t0 + 0.5 * dt);
    const auto k3 = ydot (y0 + k2 * 0.5 * dt, t0 + 0.5 * dt);
    const auto k4 = ydot (y0 + k3 * 1.0 * dt, t0 + 1.0 * dt);

    return (k1 + k2 * 2 + k3 * 2 + k4) * dt / 6;
}

std::valarray<double> RungeKuttaVector::integrateSubdivided (std::valarray<double> y0, double t0, double dt, int steps) const
{
    for (int n = 0; n < steps; ++n)
    {
        y0 += evaluateDeltaY (y0, t0, dt / steps);
        t0 += dt / steps;
    }
    return y0;
}
