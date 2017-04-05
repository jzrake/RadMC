#ifndef RandomVariable_hpp
#define RandomVariable_hpp


#include <random>
#include "TabulatedFunction.hpp"




class RandomVariable
{
public:
    class SamplingScheme
    {
    public:
        virtual double generate (double uniformSample) = 0;
        virtual ~SamplingScheme() {}
    };

    /**
    Sample the uniform random variable between 0 and 1 using the global RNG.
    */
    static double sampleUniform();

    /**
    Sample the uniform random variable between 0 and 2 pi using the global RNG.
    */
    static double sampleUniformAzimuth();

    /**
    Return a sampling scheme based on a known quantile function.
    */
    static SamplingScheme* fromQnt (std::function<double (double)> qnt);

    /**
    Return a sampling scheme from a probability density function (PDF).
    Internally, this constructs a lookup table for the cumulative distribution
    function (CDF) by numerically integrating
    */
    static SamplingScheme* fromPdf (std::function<double (double)> pdf, double x0, double x1);

    /**
    Return a "random variable" that always returns x when sampled.
    */
    static RandomVariable diracDelta (double x);

    /**
    Make a random variable which will be sampled using the given scheme.
    */
    RandomVariable (SamplingScheme* scheme);

    /**
    Make a random variable based on a known quantile function (this is just a
    shortcut for giving the return value of fromQnt to the general
    construction).
    */
    RandomVariable (std::function<double (double)> qnt);

    double sample();

private:
    // Use of shared pointer allows the random variable to be copy-constructed
    // stupidly. The underlying scheme will remain accessible to any new
    // instances.
    std::shared_ptr<SamplingScheme> scheme;
};


#endif