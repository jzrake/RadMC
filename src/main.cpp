#include <iostream>
#include "TurbulenceModelDriver.hpp"
#include "ComptonizationModelDriver.hpp"

#include <fstream>
#include "Distributions.hpp"


void testRandomVariable (RandomVariable& x)
{
    TabulatedFunction f = TabulatedFunction::makeHistogram (x.sample (1e8), 256,
        TabulatedFunction::useEqualBinWidthsLogarithmic, true, true, false);

    std::ofstream out ("x.dat");
    f.outputTable (out);
}


int main (int argc, const char *argv[])
{
    // TurbulenceModelDriver driver;
    // driver.run (argc, argv);

    ComptonizationModelDriver comptonizer;
    comptonizer.run (argc, argv);

    // double kT = 0.01;
    // auto electronPdf = Distributions::makeMaxwellJuttner (kT, Distributions::Pdf);
    // RandomVariable electronGammaBeta = RandomVariable::fromPdf (electronPdf, 0, 2);

    // testRandomVariable (electronGammaBeta);

    return 0;
}
