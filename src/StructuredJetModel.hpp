#pragma once
#include "FourVector.hpp"
#include "TabulatedFunction.hpp"
#include "RelativisticWind.hpp"




class StructuredJetModel
{
public:

    struct Config
    {
        int tableResolutionRadius = 256;
        int tableResolutionTheta = 256;
        double outermostRadius = 1e4;
        double jetOpeningAngle = 0.1;
        double jetStructureExponent = 1.0;
        double specificWindPower = 1e2;
        double luminosityPerSteradian = 1e48;
        double innerRadiusCm = 1e8;
        double leptonsPerBaryon = 1.0;
        double photonsPerBaryon = 1e4;
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
        FourVector position;
        FourVector momentum;
    };

    /** Constructor for this model. */
    StructuredJetModel (Config config);

    /**
    Get a uniformly distributed polar angle between 0 and the tabulated polar
    region.
    */
    double sampleTheta() const;

    /**
    Return an approximation of the photospheric radius (in cm) at the given
    polar angle. The estimate is very good, as long as the photosphere lies
    above the radius were the jet reaches its terminal Lorentz factor.
    */
    double approximatePhotosphere (double theta) const;

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

    /**
    Evolve a single photon from the given inner radius until it escapes to
    infinity.
    */
    std::vector<Photon> generatePhotonPath (double initialRadius, double theta);

private:
    double jetStructureEtaOfTheta (double theta) const;
    double sampleElectronGammaBeta (double kT) const;
    const TabulatedFunction& getTableForTheta (double theta) const;
    TabulatedFunction tabulateWindSolution (double rmax, double theta) const;
    void tabulateWindAllAngles (double rmax);
    RelativisticWind::WindState configureWindState (RelativisticWind::WindState state, FourVector position) const;

    Config config;
    std::vector<TabulatedFunction> tableOfSolutions;
    std::vector<double> tableOfThetas;
};
