#include "PhysicsConstants.hpp"




double PhysicsConstants::gramToErg (double massInGrams)
{
    return massInGrams * c * c;
}

double PhysicsConstants::gramToMeV (double massInGrams)
{
    const double ergsInOneMeV = 1.6021773e-6;
    return massInGrams * c * c / ergsInOneMeV;
}
