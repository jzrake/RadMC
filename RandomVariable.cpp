#include "RandomVariable.hpp"




// ========================================================================
class KnownQuantileFunction : public RandomVariable::SamplingScheme
{
public:
    KnownQuantileFunction (std::function<double (double)> quantileFunction) :
    quantileFunction (quantileFunction) {}

    double generate (double F)
    {
        return quantileFunction (F);
    }
private:
    std::function<double (double)> quantileFunction;
};




// ========================================================================
class KnownProbabilityDensityFunction : public RandomVariable::SamplingScheme
{
public:
    KnownProbabilityDensityFunction (std::function<double (double)> densityFunction, double x0, double x1)
    {
        const int numberOfTableEntries = 1024;
        const double accuracyParameter = 1e-12;
        const GaussianQuadrature gauss (8);

        tabulatedCDF = TabulatedFunction::createTabulatedIntegral (
            densityFunction, x0, x1, numberOfTableEntries,
            TabulatedFunction::useEqualBinWidthsLinear, gauss,
            accuracyParameter, true);
    }

    double generate (double F)
    {
        return tabulatedCDF.lookupArgumentValue (F);
    }
private:
    TabulatedFunction tabulatedCDF;
};




// ========================================================================
RandomVariable::SamplingScheme* RandomVariable::fromDensityFunction (
    std::function<double (double)> pdf,
    double x0, double x1)
{
    return new KnownProbabilityDensityFunction (pdf, x0, x1);
}

RandomVariable::SamplingScheme* RandomVariable::fromQuantileFunction (
    std::function<double (double)> qnt)
{
    return new KnownQuantileFunction (qnt);
}
