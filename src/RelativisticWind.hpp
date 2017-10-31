#pragma once
#include <vector>
#include "RungeKutta.hpp"
#include "FourVector.hpp"
#include "PhysicsConstants.hpp"




class RelativisticWind
{
public:
    class WindState
    {
    public:
        enum class Species { baryon, electron, photon };

        /**
        Constructor for a wind state. The first argument, reference to the
        wind model, is not kept.
        */
        WindState (const RelativisticWind& wind, double r, double u);
        WindState& setLuminosityPerSteradian (double L) { luminosityPerSteradian = L; return *this; }
        WindState& setInnerRadiusCm (double R) { innerRadiusCm = R; return *this; }
        WindState& setLeptonsPerBaryon (double Z) { leptonsPerBaryon = Z; return *this; }
        WindState& setPhotonsPerBaryon (double m) { photonsPerBaryon = m; return *this; }
        WindState& setPropagationAngle (UnitVector uhat) { propagationAngle = uhat; return *this; }

        /**
        Return the temperature kT / me c^2. This is jsut a conversion from the
        (dimensionless) specific enthalpy mu, based on the number of leptons
        and photons per bayron.
        */
        double temperature() const;

        /**
        Return the comoving number density of the given species, in 1 / cm^3.
        */
        double properNumberDensity (Species) const;

        /**
        Return the mean-free-path to free electron scattering, in cm, for a
        photon propagating in the given direction (defined as dl / dt where dt
        is the differential optical depth). This depends on the specific wind
        luminosity (eta) set by the user, and the absolute wind luminosity per
        Sr and lepton fraction set by the user. Klein-Nishina effects are
        currently ignored.
        */
        double thomsonMeanFreePath (UnitVector nhat) const;

        /**
        Return the four-velocity of the wind state. This depends on the
        propagation angle set by the user.
        */
        FourVector fourVelocity() const;

        // Wind state variables ===============================================
        double r; // radius (in units of inner boundary)
        double u; // four-velocity, u (in units of c)
        double g; // wind Lorentz factor
        double m; // specific enthalpy mu (in units of c^2, not including rest-mass)
        double p; // gas pressure (relative to density)
        double d; // gas density (equal to 1 / (r^2 u))
        /**
        Note: these variables are not all independent; they are computed when
        the data structure is initialized, so you should treat them as read-
        only variables. If you a different state then you should create a new
        one from the wind model, radius and four-velocity.
        */

    private:
        // From wind model ====================================================
        double specificWindPower;
        double initialFourVelocity;
        double adiabaticIndex;

        // Set by user ========================================================
        double luminosityPerSteradian = 1.0; // erg / s / Sr
        double innerRadiusCm = 1.0;          // inner radius (cm)
        double leptonsPerBaryon = 1.0;       // Z_{\pm}
        double photonsPerBaryon = 1.0;       // n-gamma / np
        UnitVector propagationAngle = UnitVector::zhat; // orientation of flow
        PhysicsConstants P;
    };

    RelativisticWind();

    /**
    Set the specific wind power, eta. If f is the rest-mass energy flux per
    steradian, the isotropic equivalent wind luminosity would be L = 4 pi f
    eta.
    */
    RelativisticWind& setSpecificWindPower (double eta);

    /**
    Set the wind velocity at the inner boundary.
    */
    RelativisticWind& setInitialFourVelocity (double u0);

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
    double specificWindPower = 10.0;
    double initialFourVelocity = 1.0;
    double adiabaticIndex = 4. / 3;
    mutable RungeKutta solver;
};
