#include <cmath>
#include <stdexcept>
#include "NewtonRaphesonSolver.hpp"




// ============================================================================
NewtonRaphesonSolver::NewtonRaphesonSolver (
    std::function<double (double)> F,
    std::function<double (double)> f) : F(F), f(f)
{

}

double NewtonRaphesonSolver::solve (double F0, double x0, double tolerance)
{
    Result result;
    return solve (F0, x0, tolerance, result);
}

double NewtonRaphesonSolver::solve (double F0, double x0, double tolerance, Result& result)
{
    result.numIterations = 0;
    result.succeeded = false;

    double xn = x0;

    while (true)
    {
        double Fn = F(xn);
        double fn = f(xn);

        xn = x0 - (Fn - F0) / fn;
        x0 = xn;

        if (std::fabs (F0 - Fn) < tolerance)
        {
            result.succeeded = true;
            break;
        }
        else if (++result.numIterations > 100)
        {
            throw std::runtime_error ("NewtonRaphesonSolver hit maximum iterations");
        }
    }

    return xn;
}
