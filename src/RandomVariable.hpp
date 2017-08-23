#ifndef RandomVariable_hpp
#define RandomVariable_hpp

#include <functional>
#include <vector>




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

    static double sampleUniformPitch();

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
    function (CDF) by numerically integrating the PDF. This function uses
    equal linear sampling between x0 and x1. If more flexibility is needed for
    the lookup table, then sub-class SamplingScheme, using the KnownPdf class
    found in RandomVariable.cpp as a template.
    */
    static SamplingScheme* fromPdf (std::function<double (double)> pdf, double x0, double x1);

    /**
    Return a "random variable" that always returns x when sampled.
    */
    static RandomVariable diracDelta (double x);

    /**
    Return a random variable that is uniform from x0 to x1
    */
    static RandomVariable uniformOver (double x0, double x1);

    /**
    Return an exponential random variable with time-constant beta.
    */
    static RandomVariable exponential (double beta);

    /**
    Return a power-law distributed random variable between a and b, with
    index p.
    */
    static RandomVariable powerLaw (double p, double a, double b);

    /**
    Make an empty random variable. This cannot be called, it's just here so
    there's a default constructor.
    */
    RandomVariable () {}

    /**
    Make a random variable which will be sampled using the given scheme.
    */
    RandomVariable (SamplingScheme* scheme);

    /**
    Make a random variable based on a known quantile function (this is just a
    shortcut for passing the result of fromQnt to the above constructor).
    */
    RandomVariable (std::function<double (double)> qnt);

    /**
    Generate a sample of this random variable.
    */
    double sample() const;

    /**
    Generate a bunch of samples.
    */
    std::vector<double> sample (int numberOfSamples) const;

    /**
    Print a table of values representing the PDF of the random variable. This
    is done by sampling the random variable the given number of times and
    generating a histogram.
    */
    void outputDistribution (std::ostream& stream, int numberOfSamples) const;

    /**
    Print a table of values representing the PDF of a function of the random
    variable.
    */
    void outputDistribution (std::ostream& stream, int numberOfSamples,
        std::function<double (double)> functionOfX) const;

private:
    // Use of shared pointer allows the random variable to be copy-constructed
    // stupidly. The underlying scheme will remain accessible to any new
    // instances.
    std::shared_ptr<SamplingScheme> scheme;
};


#endif