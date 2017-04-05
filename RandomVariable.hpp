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

    static SamplingScheme* fromQuantileFunction (std::function<double (double)> qnt);
    static SamplingScheme* fromDensityFunction (std::function<double (double)> pdf, double x0, double x1);

    RandomVariable (SamplingScheme* scheme) : randomVariableF (0, 1), scheme (scheme) {}

    template <class EngineType> double operator() (EngineType& engine)
    {
        double F = randomVariableF (engine);
        return scheme->generate (F);
    }

private:
    std::uniform_real_distribution<double> randomVariableF;
    std::unique_ptr<SamplingScheme> scheme;
};


#endif