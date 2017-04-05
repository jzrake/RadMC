#include <cmath>
#include <stdexcept>
#include "CubicInterpolant.hpp"




// ============================================================================
CubicInterpolant CubicInterpolant::fromValuesAndDerivatives (double xa, double xb, double fa, double fb, double ga, double gb)
{
    double c0 = (fb * std::pow (xa, 2) * (xa - 3 * xb) + xb * (gb * std::pow (xa, 2) * (-xa + xb) + xb * (fa * (3 * xa - xb) + ga * xa * (-xa + xb)))) / std::pow (xa - xb, 3);
    double c1 = (gb * std::pow (xa, 3) + xa * (-6 * fa + 6 * fb + (2 * ga + gb) * xa) * xb - (ga + 2 * gb) * xa * std::pow (xb, 2) - ga * std::pow (xb, 3)) / std::pow (xa - xb, 3);
    double c2 = (3 * fa * (xa + xb) - 3 * fb * (xa + xb) - (xa - xb) * (ga * xa + 2 * gb * xa + 2 * ga * xb + gb * xb)) / std::pow (xa - xb, 3);
    double c3 = (-2 * fa + 2 * fb + (ga + gb) * (xa - xb)) / std::pow (xa - xb, 3);
    return CubicInterpolant (c0, c1, c2, c3);
}

CubicInterpolant::CubicInterpolant() : CubicInterpolant (0, 0, 0, 0)
{

}

CubicInterpolant::CubicInterpolant (double c0, double c1, double c2, double c3) : c0 (c0), c1 (c1), c2 (c2), c3 (c3)
{

}

double CubicInterpolant::evaluate (double x)
{
    return c0 + c1 * x + c2 * x * x + c3 * x * x * x;
}
