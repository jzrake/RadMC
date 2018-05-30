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
        wind model, is not kept. The last argument, mf, is the free enthalpy
        at the wind base (as oppsed to the thermal enthalpy). By default it is
        set to zero.
        */
        WindState (double e, double u, double r, double heatingRate);
        WindState& setLuminosityPerSteradian (double L) { luminosityPerSteradian = L; return *this; }
        WindState& setInnerRadiusCm (double R) { innerRadiusCm = R; return *this; }
        WindState& setLeptonsPerBaryon (double Z) { leptonsPerBaryon = Z; return *this; }
        WindState& setPhotonsPerBaryon (double m) { photonsPerBaryon = m; return *this; }
        WindState& setPropagationAngle (UnitVector uhat) { propagationAngle = uhat; return *this; }

        /**
        Return the radius in cm.
        */
        double radius() const;

        /**
        Return the equation of motion for d-log e / d-log r.
        */
        double dLogedLogr() const;

        /**
        Return the equation of motion for d-log u / d-log r.
        */
        double dLogudLogr() const;

        /**
        Return the equation of motion for d-log w / d-log r.
        */
        double dLogwdLogr() const;

        /**
        Return the dimensionless heating rate, xi defined as dq / rho w / d-log r.
        */
        double heatingRateXi() const;

        /**
        Return the comoving number density of the given species, in 1 / cm^3.
        */
        double properNumberDensity (Species) const;

        /**
        Return the jet optical depth, ne st r / Gamma.
        */
        double jetOpticalDepth() const;

        /**
        Return the ratio of photons to protons if the photon number is given
        by thermodynamic equilibrium condition (Planck spectrum).
        */
        double blackbodyPhotonsPerProton() const;

        /**
        Return the temperature kT / me c^2. This is jsut a conversion from the
        (dimensionless) specific enthalpy mu, based on the number of leptons
        and photons per bayron.
        */
        double photonTemperature() const;

        /**
        Return the nominal Compton y-parameter, 4 tau Delta-theta.
        */
        double comptonParameter() const;

        /**
        Return Delta-theta, the amount by which the electron temperature exceeds
        the Compton temperature.
        */
        double deltaTheta() const;

        /**
        Return the electron temperature.
        */
        double electronTemperature() const;

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
        Return the comoving Thomsom scattering mean free path.
        */
        double thomsonMeanFreePathComoving() const;

        /**
        Return the nominal kinematic viscosity due to radiation: l c w / (1 + w).
        */
        double radiationViscosity() const;

        /**
        Return r / gamma.
        */
        double causallyConnectedScale() const;

        /**
        Return the four-velocity of the wind state. This depends on the
        propagation angle set by the user.
        */
        FourVector fourVelocity() const;

        // Wind state variables ===============================================
        double e; // eta, enthlapy per baryon, specficWindLuminosity
        double u; // four-velocity, u (in units of c)
        double r; // radius (in units of inner boundary)
        // Given to the constructor ^ derived here v ==========================
        double g; // wind Lorentz factor
        double h; // total specific enthalpy, 1 + m
        double w; // specific enthalpy mu_t (in units of c^2, not including rest-mass)
        double p; // gas pressure (relative to density)
        double d; // gas density (equal to 1 / (r^2 u))
        double s; // gas specific entropy
        double M2; // Mach number squared (non-relativistic)
    private:
        // Set by user ========================================================
        double luminosityPerSteradian = 1.0; // erg / s / Sr
        double innerRadiusCm = 1.0;          // inner radius (cm)
        double leptonsPerBaryon = 1.0;       // Z_{\pm}
        double photonsPerBaryon = 1.0;       // n-gamma / np
        double heatingRate = 0.0;            // delta-q / (delta-logr rho w), used in post
        UnitVector propagationAngle = UnitVector::zhat; // orientation of flow
        PhysicsConstants P;
    };

    RelativisticWind();

    /**
    Set the specific wind power, eta. If f is the rest-mass energy flux per
    steradian, the isotropic equivalent wind luminosity would be L = 4 pi f
    eta.
    */
    void setSpecificWindPower (double eta);

    /**
    Set the wind four-velocity at the inner boundary.
    */
    void setInitialFourVelocity (double u0);

    /**
    Set the heating rate (see source for definition).
    */
    void setHeatingRate (double zeta);

    /** Get the heating rate */
    double getHeatingRate() const { return heatingRate; }

    /**
    Integrate a wind profile to the given outer radius. The inner radius is
    always 1.0.
    */
    WindState integrate (double outerRadius) const;

    std::vector<WindState> integrateRange (double outerRadius) const;

    /**
    Get the wind solution at the requested radii.
    */
    std::vector<WindState> integrateTable (std::vector<double> radius) const;

private:
    void resetSolverFunction();
    double specificWindPower = 10.0; // at base only
    double initialFourVelocity = 1.0;
    double heatingRate = 0.0;
    mutable RungeKuttaVector vsolver;
};
