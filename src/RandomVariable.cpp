#include <random>
#include <cassert>
#include "RandomVariable.hpp"
#include "TabulatedFunction.hpp"




// ========================================================================
class KnownQnt : public RandomVariable::SamplingScheme
{
public:
    KnownQnt (std::function<double (double)> qnt) : qnt (qnt) {}

    double generate (double F) override
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

    double generate (double F) override
    {
        return tabulatedCDF.lookupArgumentValue (F);
    }

private:
    TabulatedFunction tabulatedCDF;
};




// ========================================================================
static std::uniform_real_distribution<double> uniform (0, 1);
static std::uniform_real_distribution<double> uniformPitch (-1, 1);
static std::uniform_real_distribution<double> uniformAzimuth (0, 2 * M_PI);
static std::mt19937 engine;

double RandomVariable::sampleUniform()
{
    return uniform (engine);
}

double RandomVariable::sampleUniformPitch()
{
    return uniformPitch (engine);
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

RandomVariable RandomVariable::uniformOver (double x0, double x1)
{
    return new KnownQnt ([=] (double F) { return x0 + (x1 - x0) * F; });
}

RandomVariable RandomVariable::exponential (double beta)
{
    return new KnownQnt ([=] (double F) { return -std::log (1 - F) * beta; });
}

RandomVariable RandomVariable::powerLaw (double p, double a, double b)
{
    return new KnownQnt ([=] (double F)
        {
            return std::pow (
                + (1 - F) * std::pow (a, 1 + p)
                + (0 + F) * std::pow (b, 1 + p),
                1 / (1 + p));
        });
}

RandomVariable::RandomVariable (SamplingScheme* scheme) : scheme (scheme)
{

}

RandomVariable::RandomVariable (std::function<double (double)> qnt) : scheme (new KnownQnt (qnt))
{

}

double RandomVariable::sample() const
{
    assert (scheme != nullptr);
    return scheme->generate (sampleUniform());
}

std::vector<double> RandomVariable::sample (int numberOfSamples) const
{
    std::vector<double> samples;

    for (int n = 0; n < numberOfSamples; ++n)
    {
        samples.push_back (sample());
    }
    return samples;
}

void RandomVariable::outputDistribution (std::ostream& stream, int numberOfSamples) const
{
    outputDistribution (stream, numberOfSamples, [] (double x) { return x; });
}

void RandomVariable::outputDistribution (std::ostream& stream, int numberOfSamples,
    std::function<double (double)> functionOfX) const
{
    std::vector<double> samples = sample (numberOfSamples);

    double xbar = 0;

    for (auto& x : samples)
    {
        x = functionOfX (x);
        xbar += x;
    }

    TabulatedFunction f = TabulatedFunction::makeHistogram (samples, 256,
        TabulatedFunction::useEqualBinWidthsLogarithmic, true, true, false);

    f.outputTable (stream);
}
