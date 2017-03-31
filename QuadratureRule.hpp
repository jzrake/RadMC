#ifndef QuadratureRule_hpp
#define QuadratureRule_hpp

#include <functional>




class QuadratureRule
{
public:
    struct EvaluationDetails
    {
        double error;
        int numberOfBinsUsed;
    };

    /**
    Estimate the integral of f(x) on the interval [a, b] using the derived
    class's quadrature rule.
    */
    virtual double integrate (std::function<double (double)> f, double a, double b) = 0;

    /**
    Computes an estimate of the definite integral of f(x) between x0 and x1
    using the given number of linearly spaced bins.
    */
    double integratePartitioned (std::function<double (double)> f, double x0, double x1, int numberOfBins);

    /**
    Computes an accurate estimate of the definite integral of f(x) between x0
    and x1. The number of subdivisions is increased until the integral value
    converges to within the accuracy parameter.
    */
    double computeDefiniteIntegral (std::function<double (double)> f, double x0, double x1, double accuracy);

    /**
    Computes an accurate estimate of the definite integral of f(x) between x0
    and x1. The number of subdivisions is increased until the integral value
    converges to within the accuracy parameter. Fills in an EvaluationDetails
    details struct so the number of iterations can be checked.
    */
    double computeDefiniteIntegral (std::function<double (double)> f, double x0, double x1, double accuracy, EvaluationDetails& details);
};




class ForwardEulerRule : public QuadratureRule
{
public:
    double integrate (std::function<double (double)> f, double a, double b) override;
};




class SimpsonRule : public QuadratureRule
{
public:
    double integrate (std::function<double (double)> f, double a, double b) override;
};




class GaussianQuadrature : public QuadratureRule
{
public:
    class Implementation;

    GaussianQuadrature (int polynomialDegree=6);
    ~GaussianQuadrature();

    void setPolynomialDegree (int degreeToUse);
    double integrate (std::function<double (double)> f, double a, double b) override;

private:
    std::unique_ptr<Implementation> implementation;
};



#endif