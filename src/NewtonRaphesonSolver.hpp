#ifndef NewtonRaphsesonSolver_hpp
#define NewtonRaphsesonSolver_hpp

#include <functional>




class NewtonRaphesonSolver
{
public:
    struct Result
    {
        int numIterations;
        bool succeeded;
    };

    NewtonRaphesonSolver (
        std::function<double (double)> F,
        std::function<double (double)> f);

    /**
    Return the argument value x at which F(x) = F0, starting with a guess value x0.
    */
    double solve (double F0, double x0, double tolerance=1e-10);

    /**
    Return the argument value x at which F(x) = F0, starting with a guess value x0.
    */
    double solve (double F0, double x0, double tolerance, Result& result);


private:
    std::function<double (double)> F; // The function itself
    std::function<double (double)> f; // The derivative of F
};

#endif
