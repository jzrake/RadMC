#include <stdexcept>
#include <cmath>
#include "RootBracketingSolver.hpp"


RootBracketingSolver::RootBracketingSolver (std::function<double (double)> f) : f(f)
{

}

double RootBracketingSolver::solve (double f0, double x0, double x1, double tolerance)
{
    Result result;
    return solve (f0, x0, x1, tolerance, result);
}

double RootBracketingSolver::solve (double f0, double x0, double x1, double tolerance, Result& result)
{
    if ((f(x0) > f0 || f(x1) < f0))
    {
        throw std::runtime_error ("The initial bracket does not enclose a root");
    }

    result.numIterations = 0;
    result.succeeded = false;
    double xmid;

    while (true)
    {
        xmid = 0.5 * (x0 + x1);

        if (std::fabs (f(xmid) - f0) < tolerance)
        {
            result.succeeded = true;
            break;
        }
        else if (result.numIterations > 100)
        {
            throw std::runtime_error ("RootBracketingSolver hit maximum iterations");
        }
        else if (f(xmid) < f0)
        {
            x0 = xmid;
        }
        else
        {
            x1 = xmid;
        }

        ++result.numIterations;
    }

    return xmid;
}
