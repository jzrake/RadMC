#ifndef TabulatedFunction_hpp
#define TabulatedFunction_hpp

#include <vector>
#include <functional>
#include "QuadratureRule.hpp"




class TabulatedFunction
{
public:
    enum BinSpacingMode
    {
        useArbitraryBinSpacing,
        useEqualBinWidthsLinear,
        useEqualBinWidthsLogarithmic,
        useEqualBinMasses,
    };

    /** Constructs an empty lookup table. */
    TabulatedFunction() {}

    /**
    Construct a lookup table for the data arrays x and y. The values of x are
    assumed to increase monotonically. If the x values are spaced linearly,
    then setting spacingMode to useEqualBinWidthsLinear will accelerate the
    lookup. If values of log(x) are spaced linearly then choose
    useEqualBinWidthsLogarithmic, and if y is spaced linearly then choose
    useEqualBinMasses.
    */
    TabulatedFunction (
        const std::vector<double>& x,
        const std::vector<double>& y,
        BinSpacingMode spacingMode=useArbitraryBinSpacing);

    /**
    Construct a lookup table containing a tabulated approximation to the
    definite integral f(x) between x0 and x1. The argument values are decided
    by the numberOfBins and binSpacingMode parameters (the table will have
    numberOfBins + 1 rows, corresponding to the bin edges). The mass in each
    bin is computed adaptively up to the accuracy parameter, using the
    QuadratureRule given. If normalized is true, then all entries are divided
    by the value in the final row of the table.
    */
    static TabulatedFunction createTabulatedIntegral (
        std::function<double(double)> f,
        double x0, double x1, int numberOfBins, BinSpacingMode spacingMode,
        const QuadratureRule& quadratureRule, double accuracy=1e-14, bool normalize=false);

    /**
    Create a tabulated function that is a histogram of samples. If density is
    true, then the table contains the mass divided by the bin width, dF / dx,
    rather than the mass dF in the bin. If normalize is true, then the total
    mass of the histogram will be \int {dF} = 1. When shift is set to true,
    the table will contain the bin minpoints and associated masses (or
    densities), and the size will be numberOfBins. Otherwise if shift is set
    to false, the table will have one more rows than the numberOfBins
    parameter, with the first column holding the location of bin edges, and
    the second column holding the bin masses to the right of the corresponding
    edge, and to the left of the one after it (hence the mass in the final bin
    will always be zero).
    */
    static TabulatedFunction makeHistogram (const std::vector<double>& samples,
        int numberOfBins, BinSpacingMode spacingMode,
        bool density=false, bool normalize=false, bool shift=false);

    /**
    Look up the function value y(x), using linear interpolation between bins.
    */
    double lookupFunctionValue (double x);

    /**
    Look up the argument value x(y), using linear interpolation between bins.
    */
    double lookupArgumentValue (double y);

    /**
    Return a lambda function that evaluates lookupFunctionValue (x). The
    closure references this, so be sure the TabulatedFunction instance
    remains alive longer than the returned lambda.
    */
    std::function<double (double)> getFunction();

    /**
    Return a lambda function that evaluates lookupArgumentValue (x). The
    closure references this, so be sure the TabulatedFunction instance
    remains alive longer than the returned lambda.
    */
    std::function<double (double)> getInverse();

    /**
    Print the table values as ASCII.
    */
    void outputTable (std::ostream& stream);

    /**
    Print the table values as ASCII, alongside a function that returns known
    exact values (helpful for testing).
    */
    void outputTable (std::ostream& stream, std::function<double (double)> exactYfunction);

private:
    std::vector<double> xdata;
    std::vector<double> ydata;
    BinSpacingMode spacingMode;
};


#endif // TabulatedFunction_hpp
