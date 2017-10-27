#pragma once
#include <functional>




class RungeKutta
{
public:
    using DerivativeFunction = std::function<double (double y, double t)>;

    RungeKutta (DerivativeFunction ydot=nullptr);

    /**
    Set the function y'(y, t) to be integrated. By default this is not.
    */
    void setFunction (DerivativeFunction ydotToUse);

    /** Set the current values of y and t. */
    void setValues (double y, double t);

    /**
    Set the error tolerance to use in choosing the number of subdivisions for
    the time interval.
    */
    void setTolerance (double relativeErrorToUse);

    /** Return the current t value */
    double getT() const;

    /** Return the current y value */
    double getY() const;

    /** Integrate up to the value t. */
    double integrate (double t);

    /**
    Return the number of divisions that were used in the time interval for the
    last call to intergate() .
    */
    int getDivisionsForLast() const;

private:
    double evaluateDeltaY (double y0, double t0, double dt) const;
    double integrateSubdivided (double y0, double t0, double dt, int steps) const;

    double yn = 0.0;
    double tn = 0.0;
    double relativeError = 1e-12;
    int divisions;
    DerivativeFunction ydot;
};
