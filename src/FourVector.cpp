#include <cmath>
#include <iostream>
#include "FourVector.hpp"
#include "RotationMatrix.hpp"
#include "LorentzBoost.hpp"



UnitVector UnitVector::xhat (0, 0);
UnitVector UnitVector::yhat (0, M_PI / 2);
UnitVector UnitVector::zhat (1, 0);

UnitVector UnitVector::sampleIsotropic()
{
    double mu = RandomVariable::sampleUniformPitch();
    double phi = RandomVariable::sampleUniformAzimuth();
    return UnitVector (mu, phi);
}

UnitVector UnitVector::normalizeFrom (double vx, double vy, double vz, bool normalized)
{
    double cosTheta = vz / (normalized ? 1 : std::sqrt (vx * vx + vy * vy + vz * vz));
    double phi = std::atan2 (vy, vx);
    return UnitVector (cosTheta, phi);
}

UnitVector::UnitVector (double pitchAngleMu, double azimuthalAnglePhi) :
pitchAngleMu (pitchAngleMu),
azimuthalAnglePhi (azimuthalAnglePhi)
{

}

void UnitVector::getCartesianComponents (double& nx, double& ny, double& nz) const
{
    const double cosTheta = pitchAngleMu;
    const double sinTheta = std::sqrt (1 - cosTheta * cosTheta);
    nx = sinTheta * std::cos (azimuthalAnglePhi);
    ny = sinTheta * std::sin (azimuthalAnglePhi);
    nz = cosTheta;
}

double UnitVector::getX() const
{
    double nx, ny, nz;
    getCartesianComponents (nx, ny, nz);
    return nx;
}

double UnitVector::getY() const
{
    double nx, ny, nz;
    getCartesianComponents (nx, ny, nz);
    return ny;
}

double UnitVector::getZ() const
{
    double nx, ny, nz;
    getCartesianComponents (nx, ny, nz);
    return nz;
}

double UnitVector::pitchAngleWith (const UnitVector& other) const
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

UnitVector UnitVector::sampleAxisymmetric (RandomVariable& pitchAngle)
{
    return sampleAxisymmetric (pitchAngle.sample());
}

UnitVector UnitVector::sampleAxisymmetric (double pitchAngle)
{
    double phi = RandomVariable::sampleUniformAzimuth();
    return UnitVector (pitchAngle, phi).withPolarAxis (*this);
}

std::ostream& operator<< (std::ostream& os, const UnitVector& nhat)
{
    double nx, ny, nz;
    nhat.getCartesianComponents (nx, ny, nz);
    os << "<unit vector> (" << nx << ", " << ny << ", " << nz << ")";
    return os;
}




// ============================================================================
FourVector::FourVector()
{
    for (int n = 0; n < 4; ++n)
    {
        components[n] = 0.0;
    }
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

FourVector FourVector::fromFourVelocity (double ux, double uy, double uz)
{
    double gm = std::sqrt (1 + (ux * ux + uy * uy + uz * uz));
    return FourVector (gm, ux, uy, uz);
}

FourVector FourVector::nullWithUnitVector (UnitVector nhat)
{
    double u[4] = {1, 0, 0, 0};
    nhat.getCartesianComponents (u[1], u[2], u[3]);
    return FourVector (u);
}

FourVector FourVector::fromGammaBetaAndUnitVector (double gammaBeta, UnitVector nhat)
{
    if (gammaBeta == 0) return FourVector (1, 0, 0, 0);
    double g = std::sqrt (1 + gammaBeta * gammaBeta);
    double nx, ny, nz;
    nhat.getCartesianComponents (nx, ny, nz);
    double u[4] = {g, gammaBeta * nx, gammaBeta * ny, gammaBeta * nz};
    return FourVector (u);
}

FourVector FourVector::fromBetaAndUnitVector (double beta, UnitVector nhat)
{
    if (beta == 0) return FourVector (1, 0, 0, 0);
    return fromGammaBetaAndUnitVector (beta / std::sqrt (1 - beta * beta), nhat);
}

FourVector FourVector::spaceLikeInDirection (double radius, UnitVector nhat)
{
    double u[4] = {0, 0, 0, 0};
    nhat.getCartesianComponents (u[1], u[2], u[3]);
    u[1] *= radius;
    u[2] *= radius;
    u[3] *= radius;
    return FourVector (u);
}

double FourVector::getTimeComponent() const
{
    return components[0];
}

const double& FourVector::operator[] (int index) const
{
    return components[index];
}

double& FourVector::operator[] (int index)
{
    return components[index];
}

double FourVector::getThreeVelocityMagnitude() const
{
    const double *u = components;
    return std::sqrt (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]) / u[0];
}

double FourVector::getThreeVelocityAlong (const UnitVector& nhat) const
{
    const double *u = components;
    double nx, ny, nz;
    nhat.getCartesianComponents (nx, ny, nz);
    return (u[1] * nx + u[2] * ny + u[3] * nz) / u[0];
}

double FourVector::radius() const
{
    const double *u = components;
    return std::sqrt (u[1] * u[1] + u[2] * u[2] + u[3] * u[3]);
}

double FourVector::theta() const
{
    return std::acos (components[3] / radius());
}

UnitVector FourVector::getUnitThreeVector() const
{
    const double *u = components;
    return UnitVector::normalizeFrom (u[1], u[2], u[3]);
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
    return -u[0] * v[0] + u[1] * v[1] + u[2] * v[2] + u[3] * v[3];
}

FourVector FourVector::operator+ (const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return FourVector (u[0] + v[0], u[1] + v[1], u[2] + v[2], u[3] + v[3]);
}

FourVector FourVector::operator- (const FourVector& other) const
{
    const double *u = components;
    const double *v = other.components;
    return FourVector (u[0] - v[0], u[1] - v[1], u[2] - v[2], u[3] - v[3]);
}

FourVector FourVector::operator* (double scalar) const
{
    const double *u = components;
    const double s = scalar;
    return FourVector (u[0] * s, u[1] * s, u[2] * s, u[3] * s);
}

FourVector FourVector::operator/ (double scalar) const
{
    const double *u = components;
    const double s = scalar;
    return FourVector (u[0] / s, u[1] / s, u[2] / s, u[3] / s);
}

FourVector& FourVector::operator+= (const FourVector& other)
{
    for (int n = 0; n < 4; ++n)
        components[n] += other.components[n];

    return *this;
}

FourVector& FourVector::operator-= (const FourVector& other)
{
    for (int n = 0; n < 4; ++n)
        components[n] -= other.components[n];

    return *this;
}

FourVector& FourVector::operator*= (double scalar)
{
    for (int n = 0; n < 4; ++n)
        components[n] *= scalar;

    return *this;
}

FourVector& FourVector::operator/= (double scalar)
{
    for (int n = 0; n < 4; ++n)
        components[n] /= scalar;

    return *this;
}

FourVector FourVector::transformedBy (const LorentzBoost& L) const
{
    return L * (*this);
}

FourVector& FourVector::transformBy (const LorentzBoost& L)
{
    return *this = transformedBy(L);
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

double FourVector::betaFromGammaBeta (double gammaBeta)
{
    return gammaBeta / std::sqrt (1 + gammaBeta * gammaBeta);
}

std::ostream& operator<< (std::ostream& os, const FourVector& u)
{
    u.printToStream (os);
    return os;
}
