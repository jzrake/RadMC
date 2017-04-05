#include "DistributionSampler.hpp"



// ============================================================================
DistributionSampler::DistributionSampler() : uniform (0, 1) {}

DistributionSampler::~DistributionSampler() {}

void DistributionSampler::computeQntFromDensity (std::function<double (double)> densityFunctionToUse, double x0, double x1)
{
    const int numberOfTableEntries = 1024;
    const double accuracyParameter = 1e-12;
    const GaussianQuadrature gauss (8);

    tabulatedCDF = TabulatedFunction::createTabulatedIntegral (
        densityFunctionToUse, x0, x1, numberOfTableEntries,
        TabulatedFunction::useEqualBinWidthsLinear, gauss,
        accuracyParameter, true);

    Qnt = tabulatedCDF.getInverse();
}

void DistributionSampler::setQnt (std::function<double (double)> QntToUse)
{
    Qnt = QntToUse;
}

std::vector<double> DistributionSampler::generateSamples (int numberOfSamples)
{
    std::vector<double> samples;

    for (int n = 0; n < numberOfSamples; ++n)
    {
        double F = uniform (engine);
        samples.push_back (Qnt (F));
    }

    return samples;
}
