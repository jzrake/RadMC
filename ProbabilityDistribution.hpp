#ifndef ProbabilityDistribution_hpp
#define ProbabilityDistribution_hpp

#include <functional>
#include <random>



class ProbabilityDistribution
{
public:
    ProbabilityDistribution (std::function<double (double)> pdf);

    double sample (const std::function<double(double)>& inverseCDF);
    
private:
    std::function<double (double)> pdf;
    std::uniform_real_distribution<double> uniform;
    std::mt19937 engine;
};


#endif
