#include "RandomVariable.hpp"




// ========================================================================
class KnownQnt : public RandomVariable::SamplingScheme
{
public:
    KnownQnt (std::function<double (double)> qnt) : qnt (qnt) {}

    double generate (double F)
    {
        return qnt (F);
    }
private:
    std::function<double (double)> qnt;
};




// ========================================================================
class KnownPdf : public RandomVariable::SamplingScheme
{
public:
    KnownPdf (std::function<double (double)> densityFunction, double x0, double x1)
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
static std::uniform_real_distribution<double> uniform (0, 1);
static std::uniform_real_distribution<double> uniformAzimuth (0, 2 * M_PI);
static std::mt19937 engine;

double RandomVariable::sampleUniform()
{
    return uniform (engine);
}

double RandomVariable::sampleUniformAzimuth()
{
    return uniformAzimuth (engine);
}

RandomVariable::SamplingScheme* RandomVariable::fromPdf (std::function<double (double)> pdf, double x0, double x1)
{
    return new KnownPdf (pdf, x0, x1);
}

RandomVariable::SamplingScheme* RandomVariable::fromQnt (std::function<double (double)> qnt)
{
    return new KnownQnt (qnt);
}

RandomVariable RandomVariable::diracDelta (double x)
{
    return new KnownQnt ([=] (double F) { return x; });
}

RandomVariable::RandomVariable (SamplingScheme* scheme) : scheme (scheme)
{

}

RandomVariable::RandomVariable (std::function<double (double)> qnt) : scheme (new KnownQnt (qnt))
{

}

double RandomVariable::sample()
{
    return scheme->generate (sampleUniform());
}
