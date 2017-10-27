#include <iostream>
#include "TurbulenceModelDriver.hpp"
#include "ComptonizationModelDriver.hpp"

#define CATCH_CONFIG_RUNNER
#define CATCH_CONFIG_FAST_COMPILE
#include "../lib/catch.hpp"




int main (int argc, const char *argv[])
{
    return Catch::Session().run (argc, argv);
}

