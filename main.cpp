#include <iostream>
#include <fstream>
#include <cmath>
#include "FourVector.hpp"
#include "LorentzBoost.hpp"
#include "QuadratureRule.hpp"
#include "ProbabilityDistribution.hpp"
#include "TabulatedFunction.hpp"
#include "CubicInterpolant.hpp"
#include "Distributions.hpp"
#include "DistributionSampler.hpp"





int testHistogram()
{
    auto pdf = Distributions::makePitchAngleGivenScattered (0.4, Distributions::probabilityDensityFunction);
    auto qnt = Distributions::makePitchAngleGivenScattered (0.4, Distributions::quantileFunction);

    DistributionSampler sampler;
    //sampler.setQuantileFunction (qnt);
    sampler.computeQuantileFunctionFromDensity (pdf, -1, 1);

    std::vector<double> samples = sampler.generateSamples (1 << 20);

    TabulatedFunction hist = TabulatedFunction::makeHistogram (samples, 128,
        TabulatedFunction::useEqualBinWidthsLinear, true, true, true);

    std::cout.flags (std::ios::fixed | std::ios::showpos);
    std::cout.precision (10);
    hist.outputTable (std::cout);

    return 0;
}


int main (int argc, char **argv)
{
    return testHistogram();
}
