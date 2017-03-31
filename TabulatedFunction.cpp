#include <cmath>
#include <iostream>
#include "TabulatedFunction.hpp"
#include "CubicInterpolant.hpp"
#include "RootBracketingSolver.hpp"



TabulatedFunction::TabulatedFunction (const std::vector<double>& x, const std::vector<double>& y, BinSpacingMode spacingMode) :
xdata (x),
ydata (y),
spacingMode (spacingMode)
{

}

TabulatedFunction TabulatedFunction::createTabulatedIntegral (std::function<double(double)> f,
    double x0, double x1, int numberOfBins, BinSpacingMode spacingMode,
    QuadratureRule& quadratureRule, double accuracy, bool normalize)
{
    const double totalMass = quadratureRule.computeDefiniteIntegral (f, x0, x1, accuracy);

    auto findNextBinLocation = [&] (double a)
    {
        switch (spacingMode)
        {
            case useArbitraryBinSpacing:
            {
                throw std::runtime_error ("Bin spacing cannot be arbitrary");
            }
            case useEqualBinWidthsLinear:
            {
                double dx = (x1 - x0) / numberOfBins;
                return a + dx;
            }
            case useEqualBinWidthsLogarithmic:
            {
                return a * std::pow (x1 / x0, 1.0 / numberOfBins);
            }
            case useEqualBinMasses:
            {
                double dF = totalMass / numberOfBins;

                auto massOfBinGivenEndpoint = [&] (double b)
                {
                    return quadratureRule.computeDefiniteIntegral (f, a, b, accuracy);
                };

                RootBracketingSolver solver (massOfBinGivenEndpoint);
                double remainingMass = quadratureRule.computeDefiniteIntegral (f, a, x1, accuracy);
                return remainingMass < dF ? x1 : solver.solve (dF, a, x1);
            }
        }
    };

    double x = x0;

    std::vector<double> argumentValue;
    std::vector<double> cumulativeMass;

    argumentValue.push_back (x);
    cumulativeMass.push_back (0);

    for (int n = 0; n < numberOfBins; ++n)
    {
        double a = x;
        double b = findNextBinLocation (a);
        double dF = quadratureRule.computeDefiniteIntegral (f, a, b, accuracy);

        argumentValue.push_back (b);
        cumulativeMass.push_back (dF + cumulativeMass[n]);

        x = b;
    }

    if (normalize)
    {
        double F1 = cumulativeMass.back();

        for (int n = 0; n < cumulativeMass.size(); ++n)
        {
            cumulativeMass[n] /= F1;
        }
    }

    return TabulatedFunction (argumentValue, cumulativeMass, spacingMode);
}

double TabulatedFunction::lookupFunctionValue (double x)
{
    auto findBinIndex = [=] ()
    {
        switch (spacingMode)
        {
            case useArbitraryBinSpacing:
            case useEqualBinMasses:
            {
                for (int n = 1; n < xdata.size(); ++n)
                {
                    if (xdata[n - 1] < x && x <= xdata[n])
                    {
                        return n;
                    }
                }
                return 0;
            }
            case useEqualBinWidthsLinear:
            {
                double x0 = xdata.front();
                double x1 = xdata.back();
                return 1 + int ((x - x0) / (x1 - x0) * (xdata.size() - 1));
            }
            case useEqualBinWidthsLogarithmic:
            {
                double L0 = std::log (xdata.front());
                double L1 = std::log (xdata.back());
                return 1 + int ((std::log (x) - L0) / (L1 - L0) * (xdata.size() - 1));
            }
        }
    };

    int n = findBinIndex();

    if (n <= 0 || xdata.size() <= n)
        throw std::runtime_error ("TabulatedFunction got out-of-range x value");

    double xa = xdata[n - 1];
    double xb = xdata[n];
    double ya = ydata[n - 1];
    double yb = ydata[n];

    return ya + (x - xa) * (yb - ya) / (xb - xa);
}

double TabulatedFunction::lookupArgumentValue (double y)
{
    auto findBinIndex = [=] ()
    {
        switch (spacingMode)
        {
            case useArbitraryBinSpacing:
            case useEqualBinWidthsLinear:
            case useEqualBinWidthsLogarithmic:
            {
                for (int n = 1; n < ydata.size(); ++n)
                {
                    if (ydata[n - 1] < y && y <= ydata[n])
                    {
                        return n;
                    }
                }
                return 0;
            }
            case useEqualBinMasses:
            {
                double y0 = ydata.front();
                double y1 = ydata.back();
                return 1 + int ((y - y0) / (y1 - y0) * (ydata.size() - 1));
            }
        }
    };

    int n = findBinIndex();

    if (n <= 0 || ydata.size() <= n)
        throw std::runtime_error ("TabulatedFunction got out-of-range y value");

    double xa = xdata[n - 1];
    double xb = xdata[n];
    double ya = ydata[n - 1];
    double yb = ydata[n];

    return xa + (y - ya) * (xb - xa) / (yb - ya);
}

std::function<double (double)> TabulatedFunction::getFunction()
{
    return [this] (double x) { return lookupFunctionValue (x); };
}

std::function<double (double)> TabulatedFunction::getInverse()
{
    return [this] (double y) { return lookupArgumentValue (y); };
}

void TabulatedFunction::outputTable (std::ostream& stream)
{
    auto dummyExactCumulativeMass = [] (double x) { return 0; };
    return outputTable (stream, dummyExactCumulativeMass);
}

void TabulatedFunction::outputTable (std::ostream& stream, std::function<double (double)> exactYfunction)
{
    for (int n = 0; n < ydata.size(); ++n)
    {
        double x = xdata[n];
        stream << n << " " << x << " " << ydata[n] << " " << exactYfunction (x) << std::endl;
    }
}

