#ifndef RootBracketingSolver_hpp
#define RootBracketingSolver_hpp



#include <functional>




class RootBracketingSolver
{
public:
    struct Result
    {
        int numIterations;
        bool succeeded;
    };

    RootBracketingSolver (std::function<double (double)> f);

    /**
    Return the argument value x at which F(x) = F0, starting with a guess interval [x0, x1].
    */
    double solve (double f0, double x0, double x1, double tolerance=1e-10);

    /**
    Return the argument value x at which F(x) = F0, starting with a guess interval [x0, x1].
    */
    double solve (double f0, double x0, double x1, double tolerance, Result& result);

private:
    std::function<double (double)> f; // The function to solve
};




#endif
