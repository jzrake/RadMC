#pragma once
#include <vector>
#include "RungeKutta.hpp"




class RelativisticWind
{
public:
    struct WindState
    {
        WindState (const RelativisticWind& wind, double r, double u);
        double r; // radius (in units of inner boundary)
        double u; // four-velocity, u (in units of c)
        double g; // wind Lorentz factor
        double m; // specific enthalpy mu (in units of c^2, not including rest-mass)
        double p; // gas pressure (relative to density)
        double d; // gas density (equal to 1 / (r^2 u))
    };

    RelativisticWind();

    /**
    Set the specific wind power, eta. If f is the rest-mass energy flux per
    steradian, the isotropic equivalent wind luminosity would be L = 4 pi f
    eta.
    */
    void setWindPower (double eta);

    /**
    Set the wind velocity at the inner boundary.
    */
    void setInitialFourVelocity (double u0);

    /**
    Integrate a wind profile to the given outer radius. The inner radius is
    always 1.0.
    */
    WindState integrate (double outerRadius) const;

    /**
    Get the wind solution at the requested radii.
    */
    std::vector<WindState> integrate (std::vector<double> radius) const;

private:
    void resetSolverFunction();
    double windPower = 10.0;
    double initialFourVelocity = 1.0;
    double adiabaticIndex = 4. / 3;
    mutable RungeKutta solver;
};
