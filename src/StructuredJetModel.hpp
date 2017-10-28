#pragma once
#include "FourVector.hpp"
#include "RelativisticWind.hpp"




class StructuredJetModel
{
public:

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
    StructuredJetModel();

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
    Electron generateElectron (const Photon& photon) const;

    /**
    Get the state of the wind at the given coordinates.
    */
    RelativisticWind::WindState sampleWind (double r, double theta) const;
    RelativisticWind::WindState sampleWind (const FourVector& position) const;

    /**
    Evolve a single photon from the given inner radius until it escapes to
    infinity.
    */
    std::vector<Photon> generatePhotonPath (double initialRadius, double theta);

private:
};
