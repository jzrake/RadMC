#ifndef RichardsonCascade_hpp
#define RichardsonCascade_hpp

#include "TabulatedFunction.hpp"



class RichardsonCascade
{
public:
    struct TimeScales
    {
        double getShortest();
        double eddyTurnoverTime;
        double viscousDampingTime;
        double comptonDragTime;
        double effectiveDampingTime;
    };

    RichardsonCascade (double kmax=1e4, int numBins=128);
    ~RichardsonCascade();

    /**
    Advance the spectral energy distribution for a time step dt.
    */
    void advance (double dt);

    /**
    Calculate relevant time scales from the spectral energy distribution.
    */
    TimeScales getTimeScales (int binIndex) const;

    /**
    Get the shortest time scale over all the length scales.
    */
    double getShortestTimeScale() const;

    /**
    Return the array of eddy turnover times at each length scale.
    */
    std::vector<double> getEddyTurnoverTime() const;

    /**
    Return the array of effective damping time scales at each length scale.
    This will be the viscous time above the photon mean free path scale,
    and the Compton time below it.
    */
    std::vector<double> getDampingTime() const;

    /**
    Return the photon mean free path scale (this is just a parameter).
    */
    double getPhotonMeanFreePathScale() const;

    /**
    Return the scale at which the eddy and viscous damping time scales would
    be equal in the absence of optically thin Compton drag.
    */
    double getFiducialViscousScale() const;

    /**
    Return the scale at which the eddy and Compton time scales would be be
    equal in the absence of optically thick photon viscosity. This scale is
    relevant when it is smaller than the photon mean free path scale.
    */
    double getFiducialComptonScale() const;

    /**
    Return the optically thin Compton power in units of the cascade power,
    assuming the fluctuating velocity at the photon mean free path scale has
    the inertial range value. When this number is small, there is no
    suppression of the cascade by Compton drag. Otherwise, the inertial range
    will be suppressed by this factor so that the Compton and eddy time scales
    are equal around the photon mean free path scale.
    */
    double getFiducialComptonPower() const;

    /** Return the total energy. */
    double getTotalEnergy() const;

    TabulatedFunction powerSpectrum;
    double cascadePower;           // should generally be 1
    double photonMeanFreePath;     // in units of the outer scale
    double radiativeEnergyDensity; // in units of rho c^2

private:
    double getEigenvalueAtEdge (int edgeIndex, double* binSpacing=nullptr) const;
    double getSignalTimeAtEdge (int edgeIndex) const;
};


#endif
