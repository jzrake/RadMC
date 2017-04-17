#include <cmath>
#include <iostream>
#include <algorithm>
#include "TabulatedFunction.hpp"
#include "CubicInterpolant.hpp"
#include "RootBracketingSolver.hpp"



TabulatedFunction::TabulatedFunction (const std::vector<double>& x, const std::vector<double>& y, BinSpacingMode spacingMode) :
xdata (x),
ydata (y),
spacingMode (spacingMode)
{

}

TabulatedFunction::TabulatedFunction (double x0, double x1, int numberOfBins, BinSpacingMode spacingMode)
{
    for (int n = 0; n < numberOfBins + 1; ++n)
    {
        switch (spacingMode)
        {
            case useEqualBinWidthsLinear:
            {
                xdata.push_back (x0 + (x1 - x0) * n / numberOfBins);
                break;
            }
            case useEqualBinWidthsLogarithmic:
            {
                xdata.push_back (x0 * std::pow (x1 / x0, double (n) / numberOfBins));
                break;
            }
            default:
            {
                throw std::runtime_error ("invalid bin spacing mode");
            }
        }
        ydata.push_back (0);
    }
}

TabulatedFunction TabulatedFunction::createTabulatedIntegral (
    std::function<double(double)> f,
    double x0, double x1, int numberOfBins, BinSpacingMode spacingMode,
    const QuadratureRule& quadratureRule, double accuracy, bool normalize)
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
                    return quadratureRule.computeDefiniteIntegral (f, a, b, accuracy); // PERF
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

        // std::cout << dF << " " << f(a) << " " << f(b) << std::endl;

        x = b;
    }

    if (normalize)
    {
        double F1 = cumulativeMass.back();

        if (F1 == 0.0)
        {
            throw std::runtime_error ("TabulatedFunction::createTabulatedIntegral: "
                "there was no mass in the CDF, consider different range");
        }

        for (int n = 0; n < cumulativeMass.size(); ++n)
        {
            cumulativeMass[n] /= F1;
        }
    }

    return TabulatedFunction (argumentValue, cumulativeMass, spacingMode);
}

TabulatedFunction TabulatedFunction::makeHistogram (
    const std::vector<double>& samples,
    int numberOfBins, BinSpacingMode spacingMode,
    bool density, bool normalize, bool shift)
{
    if (samples.size() == 0)
    {
        return TabulatedFunction (std::vector<double>(), std::vector<double>(), useArbitraryBinSpacing);
    }

    double x0 = *std::min_element (samples.begin(), samples.end()) - 1e-14; // To ensure we catch the first and last samples
    double x1 = *std::max_element (samples.begin(), samples.end()) + 1e-14;

    std::vector<double> binEdges;
    std::vector<double> binValues (numberOfBins + 1, 0.0);

    switch (spacingMode)
    {
        case useArbitraryBinSpacing:
        case useEqualBinMasses:
        {
            throw std::runtime_error ("Histogram bin spacing must be linear or logarithmic");
        }
        case useEqualBinWidthsLinear:
        {
            for (int n = 0; n < numberOfBins + 1; ++n)
            {
                binEdges.push_back (x0 + n * (x1 - x0) / numberOfBins);
            }
            break;
        }
        case useEqualBinWidthsLogarithmic:
        {
            for (int n = 0; n < numberOfBins + 1; ++n)
            {
                binEdges.push_back (x0 * std::pow (x1 / x0, double (n) / numberOfBins));
            }
            break;
        }
    }

    // Returns the index of the bin edge to the right of the sample position
    // (lowest possible value is 1).

    auto findBinIndex = [=] (double x)
    {
        switch (spacingMode)
        {
            case useEqualBinWidthsLinear:
            {
                double x0 = binEdges.front();
                double x1 = binEdges.back();
                return 0 + int ((x - x0) / (x1 - x0) * (binEdges.size() - 1));
            }
            case useEqualBinWidthsLogarithmic:
            {
                double L0 = std::log (binEdges.front());
                double L1 = std::log (binEdges.back());
                return 0 + int ((std::log (x) - L0) / (L1 - L0) * (binEdges.size() - 1));
            }
            default:
            {
                return -1;
            }
        }
    };

    double sampleMass = normalize ? 1.0 / samples.size() : 1.0;

    for (int n = 0; n < samples.size(); ++n)
    {
        int binIndex = findBinIndex (samples[n]);

        if (binIndex < 0 || binIndex >= binEdges.size() - 1)
        {
            throw std::runtime_error ("TabulatedFunction got out-of-range x value");
        }

        if (density)
        {
            binValues[binIndex] += sampleMass / (binEdges[binIndex + 1] - binEdges[binIndex]);
        }
        else
        {
            binValues[binIndex] += sampleMass;
        }
    }

    if (shift)
    {
        std::vector<double> xdata;
        std::vector<double> ydata;

        for (int n = 0; n < binEdges.size() - 1; ++n)
        {
            xdata.push_back (0.5 * (binEdges[n] + binEdges[n + 1]));
            ydata.push_back (binValues[n]);
        }

        // Even spacing would be OK to assume if spacing was
        // useEqualBinWidthsLinear, but even log spacing would not be
        // preserved because of the bin midpoint thing. Lookups on a shifted
        // histogram are unlikely to be used anyway; it's probably just going
        // to be printed or saved. So we'll play it safe and use arbitrary
        // spacing here.
        return TabulatedFunction (xdata, ydata, useArbitraryBinSpacing);
    }
    else
    {
        return TabulatedFunction (binEdges, binValues, spacingMode);
    }
}

int TabulatedFunction::size() const
{
    return xdata.size();
}

double& TabulatedFunction::operator[] (int index)
{
    if (index < 0 || index >= ydata.size())
    {
        throw std::runtime_error ("TabulatedFunction::operator[] index out of range");
    }
    return ydata[index];
}

const double& TabulatedFunction::operator[] (int index) const
{
    if (index < 0 || index >= ydata.size())
    {
        throw std::runtime_error ("TabulatedFunction::operator[] index out of range");
    }
    return ydata[index];
}

double TabulatedFunction::getBinEdge (int index) const
{
    if (index < 0 || index >= xdata.size())
    {
        throw std::runtime_error ("TabulatedFunction::getBinEdge index out of range");
    }
    return xdata[index];
}

double TabulatedFunction::getBinWidth (int index) const
{
    if (index < 0 || index >= xdata.size() - 1)
    {
        throw std::runtime_error ("TabulatedFunction::getBinWidth index out of range");
    }
    return xdata[index + 1] - xdata[index];
}

double TabulatedFunction::lookupFunctionValue (double x)
{
    auto findBinIndex = [=] () -> int
    {
        switch (spacingMode)
        {
            case useArbitraryBinSpacing:
            case useEqualBinMasses:
            {
                if (x == xdata.front())
                {
                    return 1;
                }

                auto lower = std::lower_bound (xdata.begin(), xdata.end(), x);
                return lower - xdata.begin();
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
    auto findBinIndex = [=] () -> int
    {
        switch (spacingMode)
        {
            case useArbitraryBinSpacing:
            case useEqualBinWidthsLinear:
            case useEqualBinWidthsLogarithmic:
            {
                if (y == ydata.front())
                {
                    return 1;
                }

                auto lower = std::lower_bound (ydata.begin(), ydata.end(), y);
                return lower - ydata.begin();
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
    {
        throw std::runtime_error ("TabulatedFunction got out-of-range y value");
    }

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

void TabulatedFunction::outputTable (std::ostream& stream) const
{
    for (int n = 0; n < ydata.size(); ++n)
    {
        double x = xdata[n];
        stream << n << " " << x << " " << ydata[n] << " " << std::endl;
    }
}

void TabulatedFunction::outputTable (std::ostream& stream, std::function<double (double)> exactYfunction) const
{
    for (int n = 0; n < ydata.size(); ++n)
    {
        double x = xdata[n];
        stream << n << " " << x << " " << ydata[n] << " " << exactYfunction (x) << std::endl;
    }
}

