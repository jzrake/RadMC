#include <iostream>
#include <fstream>
#include <cmath>
#include "FourVector.hpp"
#include "LorentzBoost.hpp"
#include "QuadratureRule.hpp"
#include "ProbabilityDistribution.hpp"
#include "TabulatedFunction.hpp"
#include "CubicInterpolant.hpp"







auto Normal = [] (double s)
{
    return [s] (double x)
    {
        return std::sqrt (1 / (2 * M_PI * s * s)) * std::exp (-x * x / (2 * s * s));
    };
};


auto Maxwellian = [] (double T)
{
    return [T] (double beta)
    {
        double e = beta * beta / 2 / T;
        return std::sqrt (8 / M_PI / T) * e * std::exp (-e);
    };
};


// auto MaxwellJuttner = [] (double u)
// {
//  double T = 10; // Note: for T=10, the normalization should be 1995.0396464211412
//  double g = std::sqrt (1 + u * u);
//  return u * u * std::exp (-g / T);
// };



int testHistogram()
{
    std::cout.flags (std::ios::fixed | std::ios::showpos);
    std::cout.precision (10);

    std::uniform_real_distribution<double> uniform (0, 1);
    std::mt19937 engine;
    std::vector<double> samples;

    SimpsonRule simpson;
    GaussianQuadrature gauss (8);

    TabulatedFunction tabulatedCDF1 = TabulatedFunction::createTabulatedIntegral (
        Normal (1), -5, 5, 256, TabulatedFunction::useEqualBinWidthsLinear, gauss, 1e-12, true);

    TabulatedFunction tabulatedCDF2 = TabulatedFunction::createTabulatedIntegral (
        Normal (1), -5, 5, 4096, TabulatedFunction::useEqualBinMasses, gauss, 1e-12, true);

    auto inverseCDF = tabulatedCDF2.getInverse();

    for (int n = 0; n < 100000000; ++n)
    {
        double F = uniform (engine);
        samples.push_back (inverseCDF (F));
    }

    TabulatedFunction hist = TabulatedFunction::makeHistogram (samples, 128,
        TabulatedFunction::useEqualBinWidthsLinear, true, true, true);

    hist.outputTable (std::cout);

    // tabulatedCDF2.outputTable (std::cout);

    return 0;
}

    
int main (int argc, char **argv)
{
    return testHistogram();

    SimpsonRule simpson;
    GaussianQuadrature gauss (8);

    double s0 = 3;

    if (argc > 1)
    {
        s0 = std::atof (argv[1]);
    }

    TabulatedFunction table = TabulatedFunction::createTabulatedIntegral (
        Maxwellian (1), 0, s0, 256, TabulatedFunction::useEqualBinMasses, gauss, 1e-12, true);

    std::cout.flags (std::ios::fixed | std::ios::showpos);
    std::cout.precision (10);

    // table.outputTable (std::cout);

    // for (int n = 1; n < 1024; ++n)
    // {
    //  double x = n * s0 / 1024;
    //  std::cout << x << " " << table.lookupFunctionValue (x) << std::endl;
    // }

    for (int n = 1; n < 1024; ++n)
    {
        double y = n / 1024.;
        std::cout << y << " " << table.lookupArgumentValue (y) << std::endl;
    }


    return 0;
}
