#include <cmath>
#include "FourVector.hpp"
#include "RotationMatrix.hpp"






UnitVector UnitVector::normalizeFrom (double vx, double vy, double vz)
{
    double cosTheta = vz / std::sqrt (vx * vx + vy * vy + vz * vz);
    double phi = std::atan2 (vy, vx);
    return UnitVector (cosTheta, phi);
}

void UnitVector::getCartesianComponents (double& nx, double& ny, double& nz) const
{
    const double cosTheta = pitchAngleMu;
    const double sinTheta = std::sqrt (1 - cosTheta * cosTheta);
    nx = sinTheta * std::cos (azimuthalAnglePhi);
    ny = sinTheta * std::sin (azimuthalAnglePhi);
    nz = cosTheta;
}

double UnitVector::getPitchAngleWithRespectTo (const UnitVector& other) const
{
    double u[3];
    double v[3];
    this->getCartesianComponents (u[0], u[1], u[2]);
    other.getCartesianComponents (v[0], v[1], v[2]);
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

UnitVector UnitVector::withPolarAxis (const UnitVector& newPolarAxis)
{
    double theta = std::acos (newPolarAxis.pitchAngleMu);
    double phi = newPolarAxis.azimuthalAnglePhi;
    return RotationMatrix::aboutZ (phi) * (RotationMatrix::aboutY (theta) * (*this));
}

std::ostream& operator<< (std::ostream& os, const UnitVector& nhat)
{
    double nx, ny, nz;
    nhat.getCartesianComponents (nx, ny, nz);
    os << "<unit vector> (" << nx << ", " << ny << ", " << nz << ")";
    return os;
}


std::uniform_real_distribution<double> UnitVector::randomVariableMu (-1, 1);
std::uniform_real_distribution<double> UnitVector::randomVariablePhi (0, 2 * M_PI);



FourVector::FourVector()
{

}

FourVector::FourVector (const FourVector& other)
{
    for (int n = 0; n < 4; ++n)
    {
        components[n] = other.components[n];
    }
}

FourVector::FourVector (double u[4])
{
    for (int n = 0; n < 4; ++n)
    {
        components[n] = u[n];
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

FourVector FourVector::nullWithUnitVector (UnitVector nhat)
{
    double u[4] = {1, 0, 0, 0};
    nhat.getCartesianComponents (u[1], u[2], u[3]);
    return FourVector (u);
}

FourVector FourVector::fromGammaBetaAndUnitVector (double gammaBeta, UnitVector nhat)
{
    double g = std::sqrt (1 + gammaBeta * gammaBeta);
    double nx, ny, nz;
    nhat.getCartesianComponents (nx, ny, nz);
    double u[4] = {g, gammaBeta * nx, gammaBeta * ny, gammaBeta * nz};
    return FourVector (u);
}

FourVector FourVector::fromBetaAndUnitVector (double beta, UnitVector nhat)
{
    return fromGammaBetaAndUnitVector (beta / std::sqrt (1 - beta * beta), nhat);
}

double FourVector::getTimeComponent() const
{
    return components[0];
}

double FourVector::getThreeVelocityMagnitude() const
{
    const double *u = components;
    return std::sqrt (1 - 1 / (u[0] * u[0]));
}

UnitVector FourVector::getUnitThreeVector() const
{
    const double *u = components;
    return UnitVector::normalizeFrom (u[1], u[2], u[3]);
}

FourVector FourVector::operator+(const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return FourVector (u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}

FourVector FourVector::operator-(const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return FourVector (u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}

FourVector FourVector::operator-() const
{
    const double *u = components;
    return FourVector (u[0], -u[1], -u[2], -u[3]);
}

FourVector FourVector::operator* (double scalar) const
{
    const double *u = components;
    const double s = scalar;
    return FourVector (u[0] * s, u[1] * s, u[2] * s, u[3] * s);
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

void FourVector::printToStream (std::ostream& stream) const
{
    stream << "<four vector> ("
    << components[0] << ", "
    << components[1] << ", "
    << components[2] << ", "
    << components[3] << ")";
}

std::ostream& operator<< (std::ostream& os, const FourVector& u)
{
    u.printToStream (os);
    return os;
}
