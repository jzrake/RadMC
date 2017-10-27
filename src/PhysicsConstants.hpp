#pragma once



/**
A class that holds common physics constants in CGS units and does mundane unit
conversions.
*/
class PhysicsConstants
{
public:
    const double c = 2.99792458e+10; // light speed
    const double k = 1.38065040e-16; // Boltzmann's constant
    const double h = 6.62606885e-27; // Planck constant (h, not h-bar)
    const double G = 6.67408310e-08; // Newton gravitational constant
    const double e = 4.80320427e-10; // electron charge (esu)
    const double me = 9.10938215e-28; // electron mass (g)
    const double mp = 1.67262190e-24; // proton mass (g)
    const double st = 6.65245871e-25; // Thomson cross section (cm^2)

    /** Convert a mass in grams to m c^2 in erg. */
    double gramToErg (double massInGrams);
    /** Convert a mass in grams to m c^2 in MeV. */
    double gramToMeV (double massInGrams);
};
