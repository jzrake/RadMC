
#ifndef FourVector_hpp
#define FourVector_hpp


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

    /**
    Return the Lorentz factor assuming this is space-like vector.
    */
    double getLorentzFactor() const;

    /**
    Return the magnitude of the vector's associated three-velocity, assuming
    this is a space-like vector.
    */
    double getThreeVelocityMagnitude() const;

    /**
    Return four vector whose spatial components are negated.
    */
    FourVector operator-() const;

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

private:
    friend class LorentzBoost;
    double components[4];
};


#endif
