#include "ProbabilityDistribution.hpp"
#include "NewtonRaphesonSolver.hpp"
#include "QuadratureRule.hpp"



ProbabilityDistribution::ProbabilityDistribution (std::function<double (double)> pdf) :
pdf (pdf),
uniform (0, 1)
{

}

double ProbabilityDistribution::sample (const std::function<double(double)>& inverseCDF)
{
    return inverseCDF (uniform (engine));
}
