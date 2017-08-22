#include <cmath>
#include <stdexcept>
#include <iostream>
#include "QuadratureRule.hpp"




// ============================================================================
double QuadratureRule::integratePartitioned (std::function<double (double)> f, double x0, double x1,
    int numberOfBins) const
{
    double dx = (x1 - x0) / numberOfBins;
    double x = x0;
    double F = 0.0;

    for (int n = 0; n < numberOfBins; ++n)
    {
        double a = x;
        double b = x + dx;
        double dF = integrate (f, a, b);

        x += dx;
        F += dF;
    }

    return F;
}

double QuadratureRule::computeDefiniteIntegral (std::function<double (double)> f, double x0, double x1,
    double accuracy) const
{
    EvaluationDetails details;
    return computeDefiniteIntegral (f, x0, x1, accuracy, details);
}

double QuadratureRule::computeDefiniteIntegral (std::function<double (double)> f, double x0, double x1,
    double accuracy, EvaluationDetails& details) const
{
    int N = 1;
    double F0;
    double F1;
    double error;

    do
    {
        F0 = integratePartitioned (f, x0, x1, N);
        F1 = integratePartitioned (f, x0, x1, N * 2);

        error = std::fabs ((F1 - F0) / F1);
        N *= 2;

        if (N > 1<<14)
        {
            throw std::runtime_error ("QuadratureRule not making progress toward accuracy goal");
        }

        // std::cout << N << " " << f(x0) << " " << f(x1) << " " << error << std::endl;

    } while (error > accuracy);

    details.error = error;
    details.numberOfBinsUsed = N;

    return F1;
}




// ============================================================================
double ForwardEulerRule::integrate (std::function<double (double)> f, double a, double b)const
{
    return (b - a) * f(a);
}




// ============================================================================
double SimpsonRule::integrate (std::function<double (double)> f, double a, double b)const
{
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b));
}




// ============================================================================
class GaussianQuadrature::Implementation
{
public:
    virtual double integrate (std::function<double (double)> f, double a, double b) const = 0;
};




// https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
// ============================================================================
template <int N> class GaussLegendreQuadrature : public GaussianQuadrature::Implementation
{
public:
    enum {eDEGREE = N};

    double integrate (std::function<double (double)> f, double a, double b) const override
    {
        const LegendrePolynomial& legpoly = legendrePolynomial;

        double p = (b - a) / 2;
        double q = (b + a) / 2;
        double sum = 0;

        for (int i = 1; i <= eDEGREE; ++i)
        {
            sum += legpoly.weight (i) * f (p * legpoly.root(i) + q);
        }

        return p * sum;
    }

private:

    class LegendrePolynomial
    {
    public:
        LegendrePolynomial()
        {
            for (int i = 0; i <= eDEGREE; ++i)
            {
                double dr = 1;

                // Find zero
                Evaluation eval (std::cos (M_PI * (i - 0.25) / (eDEGREE + 0.5)));
                do
                {
                    dr = eval.v() / eval.d();
                    eval.evaluate (eval.x() - dr);
                } while (std::fabs (dr) > 2e-16);

                this->_r[i] = eval.x();
                this->_w[i] = 2 / ((1 - eval.x() * eval.x()) * eval.d() * eval.d());
            }
        }

        double root (int i) const { return this->_r[i]; }
        double weight (int i) const { return this->_w[i]; }

    private:
        double _r[eDEGREE + 1];
        double _w[eDEGREE + 1];

        class Evaluation
        {
        public:
            explicit Evaluation (double x) : _x(x), _v(1), _d(0)
            {
                this->evaluate (x);
            }

            void evaluate (double x)
            {
                this->_x = x;

                double vsub1 = x;
                double vsub2 = 1;
                double f     = 1 / (x * x - 1);

                for (int i = 2; i <= eDEGREE; ++i) {
                    this->_v = ((2 * i - 1) * x * vsub1 - (i - 1) * vsub2) / i;
                    this->_d = i * f * (x * this->_v - vsub1);

                    vsub2 = vsub1;
                    vsub1 = this->_v;
                }
            }

            double v() const { return this->_v; }
            double d() const { return this->_d; }
            double x() const { return this->_x; }

        private:
            double _x;
            double _v;
            double _d;
        };
    };

    static LegendrePolynomial legendrePolynomial; // Pre-compute the weights and abscissae of the Legendre polynomials
};

template <int N> typename GaussLegendreQuadrature<N>::LegendrePolynomial GaussLegendreQuadrature<N>::legendrePolynomial;




// ============================================================================
GaussianQuadrature::GaussianQuadrature (int polynomialDegree)
{
    setPolynomialDegree (polynomialDegree);
}

GaussianQuadrature::~GaussianQuadrature()
{

}

void GaussianQuadrature::setPolynomialDegree (int degreeToUse)
{
    switch (degreeToUse)
    {
        case 2: implementation.reset (new GaussLegendreQuadrature<2>()); break;
        case 3: implementation.reset (new GaussLegendreQuadrature<3>()); break;
        case 4: implementation.reset (new GaussLegendreQuadrature<4>()); break;
        case 5: implementation.reset (new GaussLegendreQuadrature<5>()); break;
        case 6: implementation.reset (new GaussLegendreQuadrature<6>()); break;
        case 7: implementation.reset (new GaussLegendreQuadrature<7>()); break;
        case 8: implementation.reset (new GaussLegendreQuadrature<8>()); break;
        case 9: implementation.reset (new GaussLegendreQuadrature<9>()); break;
        case 10: implementation.reset (new GaussLegendreQuadrature<10>()); break;
        case 11: implementation.reset (new GaussLegendreQuadrature<11>()); break;
        case 12: implementation.reset (new GaussLegendreQuadrature<12>()); break;
        case 13: implementation.reset (new GaussLegendreQuadrature<13>()); break;
        case 14: implementation.reset (new GaussLegendreQuadrature<14>()); break;
        case 15: implementation.reset (new GaussLegendreQuadrature<15>()); break;
        case 16: implementation.reset (new GaussLegendreQuadrature<16>()); break;
        case 17: implementation.reset (new GaussLegendreQuadrature<17>()); break;
        case 18: implementation.reset (new GaussLegendreQuadrature<18>()); break;
        case 19: implementation.reset (new GaussLegendreQuadrature<19>()); break;
        case 20: implementation.reset (new GaussLegendreQuadrature<20>()); break;
        case 21: implementation.reset (new GaussLegendreQuadrature<21>()); break;
        case 22: implementation.reset (new GaussLegendreQuadrature<22>()); break;
        case 23: implementation.reset (new GaussLegendreQuadrature<23>()); break;
        case 24: implementation.reset (new GaussLegendreQuadrature<24>()); break;
    }
}

double GaussianQuadrature::integrate (std::function<double (double)> f, double a, double b) const
{
    return implementation->integrate (f, a, b);
}
