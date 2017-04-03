#include <cmath>
#include <stdexcept>
#include "Distributions.hpp"




std::function<double (double)> Distributions::makeGaussian (double meanMu, double stdDeviationSigma, FunctionType type)
{
    switch (type)
    {
        case probabilityDensityFunction:
        {
            return [=] (double x)
            {
                const double m = meanMu;
                const double s = stdDeviationSigma;
                return std::sqrt (1 / (2 * M_PI * s * s)) * std::exp (-(x - m) * (x - m) / (2 * s * s));
            };
        }
        case probabilityMassFunction:
        {
            return [=] (double x)
            {
                const double m = meanMu;
                const double s = stdDeviationSigma;
                return 0.5 * (1 + std::erf ((x - m) / (s * std::sqrt (2))));
            };
        }
        case quantileFunction:
        {
            throw std::runtime_error ("not implemented");
        }
    }
}

std::function<double (double)> Distributions::makeMaxwellBoltzmann (double temperatureTheta, FunctionType type)
{
    switch (type)
    {
        case probabilityDensityFunction:
        {
            return [=] (double beta)
            {
                const double T = temperatureTheta;
                const double e = beta * beta / 2 / T;
                return std::sqrt (8 / M_PI / T) * e * std::exp (-e);
            };
        }
        case probabilityMassFunction:
        case quantileFunction:
        {
            throw std::runtime_error ("not implemented");
        }
    }
}

std::function<double (double)> Distributions::makeMaxwellJuttner (double temperatureTheta, FunctionType type)
{
    // Note: for T=10, the normalization should be 1995.0396464211412
    switch (type)
    {
        case probabilityDensityFunction:
        {
            return [=] (double gammaBeta)
            {
                const double T = temperatureTheta;
                const double u = gammaBeta;
                const double g = std::sqrt (1 + u * u);
                return u * u * std::exp (-g / T);
            };
        }
        case probabilityMassFunction:
        case quantileFunction:
        {
            throw std::runtime_error ("not implemented");
        }
    }
}

std::function<double (double)> Distributions::makePitchAngleGivenScattered (double velocityBeta, FunctionType type)
{
    switch (type)
    {
        case probabilityDensityFunction:
        {
            return [=] (double pitchAngleMu)
            {
                return 0.5 * (1 - velocityBeta * pitchAngleMu);
            };
        }
        case probabilityMassFunction:
        {
            return [=] (double pitchAngleMu)
            {
                return 0.5 * (1 - velocityBeta * pitchAngleMu);
            };
        }
        case quantileFunction:
        {
            return [=] (double quantileF)
            {
                const double B = velocityBeta;
                const double a = -B / 2;
                const double b = 1;
                const double c = 1 + B / 2 - 2 * quantileF;
                const double mu = (-b + std::sqrt (b * b - 4 * a * c)) / (2 * a);
                return mu;
            };
        }
    }
}

