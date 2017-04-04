#ifndef RandomVariable_hpp
#define RandomVariable_hpp


#include <random>
#include "TabulatedFunction.hpp"


class RandomVariable
{
public:

    // ========================================================================
    class FromQuantileFunction
    {
    public:
        FromQuantileFunction (std::function<double (double)> quantileFunction) :
        randomVariableF (0, 1),
        quantileFunction (quantileFunction)
        {

        }
        
        template <class EngineType> double operator() (EngineType& engine)
        {
            return quantileFunction (randomVariableF (engine));
        }

    private:
        std::uniform_real_distribution<double> randomVariableF;
        std::function<double (double)> quantileFunction;
    };

    // ========================================================================
    class FromProbabilityDensityFunction
    {
    public:
        FromProbabilityDensityFunction (std::function<double (double)> densityFunctionToUse, double x0, double x1) :
        randomVariableF (0, 1)
        {
            const int numberOfTableEntries = 1024;
            const double accuracyParameter = 1e-12;
            const GaussianQuadrature gauss (8);

            tabulatedCDF = TabulatedFunction::createTabulatedIntegral (
                densityFunctionToUse, x0, x1, numberOfTableEntries,
                TabulatedFunction::useEqualBinWidthsLinear, gauss,
                accuracyParameter, true);
        }

        template <class EngineType> double operator() (EngineType& engine)
        {
            return tabulatedCDF.lookupArgumentValue (randomVariableF (engine));
        }

    private:
        std::uniform_real_distribution<double> randomVariableF;
        TabulatedFunction tabulatedCDF;
    };
};


#endif