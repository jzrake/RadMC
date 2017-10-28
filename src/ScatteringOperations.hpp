#pragma once
#include "FourVector.hpp"




class ScatteringOperations
{
public:
    ScatteringOperations();

    /**
    Sample the probability space of particles that scatter with a photon whose
    propagation vector is k. The particles are mono-energetic (four velocity
    u) and isotropic in the system frame.
    */
    FourVector sampleScatteredParticles (FourVector k, double u) const;

    /**
    Sample the probability space of particles that scatter with a photon whose
    propagation vector is k. The particles are mono-energetic (four velocity
    u) and isotropic in the frame moving with the four velocity indicated by
    the restFrame argument.
    */
    FourVector sampleScatteredParticles (FourVector restFrame, FourVector k, double u) const;

    /**
    Perform a Compton scattering operation between a photon and electron with
    the given momenum four-vectors. Return the momentum change of the photon.
    */
    FourVector comptonScatter (FourVector& photon, FourVector& electron) const;


private:
};
