#ifndef DistributionSampler_hpp
#define DistributionSampler_hpp

#include <vector>
#include <random>
#include "TabulatedFunction.hpp"



class DistributionSampler
{
public:
    DistributionSampler();
    virtual ~DistributionSampler();
    
    virtual void computeQntFromDensity (std::function<double (double)> densityFunctionToUse, double x0, double x1);
    void setQnt (std::function<double (double)> QntToUse);
    std::vector<double> generateSamples (int numberOfSamples);

private:
    std::uniform_real_distribution<double> uniform;
    std::mt19937 engine;
    std::function<double (double)> Qnt;
    TabulatedFunction tabulatedCDF;
};


#endif
