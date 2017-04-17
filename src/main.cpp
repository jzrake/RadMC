#include <iostream>
#include "TurbulenceModelDriver.hpp"
#include "ComptonizationModelDriver.hpp"



int main (int argc, const char *argv[])
{
    // TurbulenceModelDriver driver;
    // driver.run (argc, argv);


    ComptonizationModelDriver comptonizer;
    comptonizer.run (argc, argv);


    return 0;
}
