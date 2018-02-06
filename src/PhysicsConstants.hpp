#pragma once



/**
A class that holds common physics constants in CGS units and does mundane unit
conversions.
*/
class PhysicsConstants
{
public:
    double c = 2.99792458e+10; // light speed
    double k = 1.38065040e-16; // Boltzmann's constant
    double h = 6.62606885e-27; // Planck constant (h, not h-bar)
    double G = 6.67408310e-08; // Newton gravitational constant
    double e = 4.80320427e-10; // electron charge (esu)
    double me = 9.10938215e-28; // electron mass (g)
    double mp = 1.67262190e-24; // proton mass (g)
    double st = 6.65245871e-25; // Thomson cross section (cm^2)
    double pc = 3.08567757e+18; // parsec (cm)
    double pi = 3.14159265e+00; // pi
    /** Convert a mass in grams to m c^2 in erg. */
    double gramToErg (double massInGrams);
    /** Convert a mass in grams to m c^2 in MeV. */
    double gramToMeV (double massInGrams);
};
