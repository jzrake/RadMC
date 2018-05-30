#pragma once
#include "FourVector.hpp"
#include "TabulatedFunction.hpp"
#include "RelativisticWind.hpp"
#include "PhysicsConstants.hpp"




class StructuredJetModel
{
public:

    struct Config
    {
        bool approximateElectronEnergiesAsDelta = true;
        int tableResolutionRadius = 256;
        int tableResolutionTheta = 256;
        double outermostRadius = 1e4;
        double jetOpeningAngle = 0.1;
        double jetPolarBoundary = 0.5;
        double jetStructureExponent = 1.0;
        double specificWindPower = 1e2;
        double luminosityPerSteradian = 1e48;
        double heatingRate = 0.0;
        double innerRadiusCm = 1e8;
        double leptonsPerBaryon = 1.0;
        double photonsPerBaryon = 1e4;
    };

    enum class Variable
    {
        FourVelocity,
        SpecificWindPower,
    };

    /**
    The electron data structure for this model.
    */
    class Electron
    {
    public:
        Electron () : momentum (1, 0, 0, 0) {}
        Electron (FourVector momentum) : momentum (momentum) {}
        FourVector momentum;
    };

    /**
    The photon data structure for this model.
    */
    class Photon
    {
    public:
        Photon (FourVector momentum=FourVector()) : momentum (momentum) {}

        /**
        Return the difference in time (seconds) between the moment a distant
        observer receives this photon, and the time when he would receive a
        light pulse which left the origin at t=0. The observer is located
        along the photon's propagation vector, relative to the origin, not the
        photon position vector. The parameter lengthUnits is the radius in cm
        when position.radius() == 1.
        */
        double lagTime (double lengthUnits) const;

        FourVector position;
        FourVector momentum;
    };

    /** Constructor for this model. */
    StructuredJetModel (Config config);

    /**
    Return a theta value, sampled uniformly over the given fraction of solid
    angle covered by the tabulated solution. Setting the fraction too close to
    1.0 may result in many photons diffusing out the sides of the jet.
    */
    double sampleTheta (double fraction) const;

    /**
    Return an approximation of the photospheric radius (in cm) at the given
    polar angle. The estimate is very good, as long as the photosphere lies
    above the radius were the jet reaches its terminal Lorentz factor.
    */
    double approximatePhotosphere (double theta) const;

    /**
    Return an approximate lag time (in seconds) for emission into the given
    polar angle.
    */
    double approximateLagTime (double theta) const;

    /**
    Return the luminosity (erg/s/Sr) at the give polar angle.
    */
    double angularLuminosity (double theta) const;

    /**
    Return the total luminosity (in erg/s) for this jet structure.
    */
    double totalLuminosity() const;

    /**
    Return the time (in seconds, source frame) for a fluid element to propagte
    from the wind inner boundary to the given radius. NOTE: the input radius
    is in units of the inner boundary, not in cm.
    */
    double fluidPropagationTimeToRadius (double r, double theta) const;

    /**
    Advance the given photon by first scattering it, and then moving it by a
    single scattering length (with respect to the initial position).
    */
    Photon stepPhoton (const Photon& photon) const;

    /**
    Generate a new photon at the given radius and polar angle. The photon is
    sampled at a random azimuth phi, and isotropically in the wind rest frame
    at that point. The photon energy is exactly the local temperature.
    */
    Photon generatePhoton (double r, double theta) const;

    /**
    Generate a new electron to scatter with the given photon based on its
    location and propagation direction. The electron energy is sampled from a
    relativistic Maxwellian with the local wind temperature.
    */
    Electron generateElectron (const Photon& photon, RelativisticWind::WindState state) const;

    /**
    Get the state of the wind at the given coordinates.
    */
    RelativisticWind::WindState sampleWindSpherical (double r, double theta) const;

    /**
    Get the state of the wind at the given coordinates.
    */
    RelativisticWind::WindState sampleWind (const FourVector& position) const;

private:
    double jetStructureEtaOfTheta (double theta) const;
    double jetStructureEffOfTheta (double theta) const;
    double sampleElectronGammaBeta (double kT) const;
    const TabulatedFunction& getTableForTheta (double theta, Variable var) const;
    TabulatedFunction tabulateWindSolution (double rmax, double theta, Variable var) const;
    void tabulateWindAllAngles (double rmax);
    RelativisticWind::WindState configureWindState (RelativisticWind::WindState state, FourVector position) const;
    RelativisticWind makeWindSolver (double theta) const;

    Config config;
    std::vector<TabulatedFunction> tableOfSolutionsU;
    std::vector<TabulatedFunction> tableOfSolutionsE;
    std::vector<double> tableOfThetas;
    PhysicsConstants physics;
};
