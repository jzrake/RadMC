#include <cmath>
#include "FourVector.hpp"


FourVector::FourVector()
{

}

FourVector::FourVector (const FourVector& other)
{
    for (int i = 0; i < 4; ++i)
    {
        components[i] = other.components[i];
    }
}

FourVector::FourVector (double u[4])
{
    for (int i = 0; i < 4; ++i)
    {
        components[i] = u[i];
    }
}

FourVector::FourVector (double E, double px, double py, double pz)
{
    components[0] = E;
    components[1] = px;
    components[2] = py;
    components[3] = pz;
}

FourVector FourVector::fromThreeVelocity (double vx, double vy, double vz)
{
    double gm = 1.0 / std::sqrt (1 - (vx * vx + vy * vy + vz * vz));
    return FourVector (gm, gm * vx, gm * vy, gm * vz);
}

double FourVector::getLorentzFactor() const
{
    return components[0];
}

double FourVector::getThreeVelocityMagnitude() const
{
    return sqrt (1 - 1 / (components[0] * components[0]));
}

FourVector FourVector::operator-() const
{
    const double *u = components;
    return FourVector (u[0], -u[1], -u[2], -u[3]);
}

double FourVector::operator* (const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return -v[0] * v[0] + u[1] * v[1] + u[2] * v[2] + u[3] * v[3];
}

bool FourVector::isNull (double tol) const
{
    return std::fabs ((*this) * (*this)) < tol;
}

bool FourVector::isFourVelocity (double tol) const
{
    return std::fabs ((*this) * (*this) + 1) < tol;
}

bool FourVector::isSpacelike() const
{
    return (*this) * (*this) < 0;
}

bool FourVector::isTimelike() const
{
    return (*this) * (*this) > 0;
}
