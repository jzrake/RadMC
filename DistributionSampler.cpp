#include "DistributionSampler.hpp"



// ============================================================================
DistributionSampler::DistributionSampler() : uniform (0, 1) {}

DistributionSampler::~DistributionSampler() {}

void DistributionSampler::computeQuantileFunctionFromDensity (std::function<double (double)> densityFunctionToUse, double x0, double x1)
{
    const int numberOfTableEntries = 1024;
    const double accuracyParameter = 1e-12;
    const GaussianQuadrature gauss (8);

    tabulatedCDF = TabulatedFunction::createTabulatedIntegral (
        densityFunctionToUse, x0, x1, numberOfTableEntries,
        TabulatedFunction::useEqualBinWidthsLinear, gauss,
        accuracyParameter, true);

    quantileFunction = tabulatedCDF.getInverse();
}

void DistributionSampler::setQuantileFunction (std::function<double (double)> quantileFunctionToUse)
{
    quantileFunction = quantileFunctionToUse;
}

std::vector<double> DistributionSampler::generateSamples (int numberOfSamples)
{
    std::vector<double> samples;

    for (int n = 0; n < numberOfSamples; ++n)
    {
        double F = uniform (engine);
        samples.push_back (quantileFunction (F));
    }

    return samples;
}
