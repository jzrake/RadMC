#ifndef Distributions_hpp
#define Distributions_hpp

#include <functional>

class Distributions
{
public:
    enum FunctionType
    {
        probabilityDensityFunction = 1,
        probabilityMassFunction    = 2,
        quantileFunction           = 3,
    };
    static std::function<double (double)> makeGaussian (double meanMu, double stdDeviationSigma, FunctionType type);
    static std::function<double (double)> makeMaxwellBoltzmann (double temperatureTheta, FunctionType type);
    static std::function<double (double)> makeMaxwellJuttner (double temperatureTheta, FunctionType type);
    static std::function<double (double)> makePitchAngleGivenScattered (double velocityBeta, FunctionType type);
};

#endif