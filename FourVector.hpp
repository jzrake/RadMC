
#ifndef FourVector_hpp
#define FourVector_hpp

#include <ostream>
#include <random>



/**
A class to encapsulate a 3D unit vector.
*/
class UnitVector
{
public:
    template <class EngineType> static UnitVector generateIsotropic (EngineType& engine)
    {
        return UnitVector (randomVariableMu (engine), randomVariablePhi (engine));
    }
    static UnitVector normalizeFrom (double vx, double vy, double vz);

    UnitVector (double pitchAngleMu, double azimuthalAnglePhi) :
    pitchAngleMu (pitchAngleMu),
    azimuthalAnglePhi (azimuthalAnglePhi) {}

    /**
    Return the cartesian components of this unit vector.
    */
    void getCartesianComponents (double& nx, double& ny, double& nz) const;

    /**
    Return the cosine of the angle between two unit vectors.
    */
    double getPitchAngleWithRespectTo (const UnitVector& other) const;

    /**
    Return the unit vector which results if this vector's pitch and azimuthal
    angles were with respect to newPolarAxis, rather than the z-axis. For
    example: zhat.withPolarAxis (xhat) = xhat.
    */
    UnitVector withPolarAxis (const UnitVector& newPolarAxis);

    double pitchAngleMu;
    double azimuthalAnglePhi;

private:
    static std::uniform_real_distribution<double> randomVariableMu;
    static std::uniform_real_distribution<double> randomVariablePhi;
};


std::ostream& operator<< (std::ostream& os, const UnitVector& nhat);




/**
A class to encapsulate four-vector operations. The (-,+,+,+) metric is assumed.
*/
class FourVector
{
public:
    FourVector();
    FourVector (const FourVector& other);
    FourVector (double u[4]);
    FourVector (double E, double px, double py, double pz);
    static FourVector fromThreeVelocity (double vx, double vy, double vz);
    static FourVector nullWithUnitVector (UnitVector nhat);
    static FourVector fromGammaBetaAndUnitVector (double gammaBeta, UnitVector nhat);
    static FourVector fromBetaAndUnitVector (double beta, UnitVector nhat);

    /**
    Return the magnitude of the vector's associated three-velocity, sqrt (1 /
    u0^2), where u0 is the Lorentz factor. This assumes isFourVelocity() would
    evaluate to true.
    */
    double getThreeVelocityMagnitude() const;

    /**
    Return a unit vector from this four vector's spatial components.
    */
    UnitVector getUnitThreeVector() const;

    /**
    Return the time component u[0] of the four vector.
    */
    double getTimeComponent() const;

    /**
    Return four vector addition.
    */
    FourVector operator+ (const FourVector& other) const;

    /**
    Return four vector subtraction.
    */
    FourVector operator- (const FourVector& other) const;

    /**
    Return four vector whose spatial components are negated.
    */
    FourVector operator-() const;

    /**
    Return four vector with all components multiplied by the given scalar.
    */
    FourVector operator* (double scalar) const;

    /**
    Return the contraction of this with another four-vector u.
    */
    double operator* (const FourVector& other) const;

    /**
    Return true if this four vector has zero contraction with itself.
    */
    bool isNull (double tol=1e-10) const;

    /**
    Return true if this four vector is a four-velocity, that is u.u == -1.
    */
    bool isFourVelocity (double tol=1e-10) const;

    /**
    Return true if this four vector has positive contraction with itself.
    */
    bool isSpacelike() const;

    /**
    Return true if this four vector has negative contraction with itself.
    */
    bool isTimelike() const;

    /**
    Print the four vector in a format like <four vector>: (ut, ux, uy, uz).
    */
    void printToStream (std::ostream& stream) const;

private:
    friend class LorentzBoost;
    double components[4];
};

std::ostream& operator<< (std::ostream& os, const FourVector& u);


#endif
