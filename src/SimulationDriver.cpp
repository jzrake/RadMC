#include <iostream>
#include <iomanip>
#include "SimulationDriver.hpp"
#include "PathHelpers.hpp"




// ============================================================================
SimulationDriver::Status::Status()
{
    outputsWrittenSoFar = 0;
    simulationIter = 0;
    simulationTime = 0.0;
}




// ============================================================================
void SimulationDriver::run (int argc, const char *argv[])
{
    makeUserParameters (userParams);

    try
    {
        Variant::updateFromCommandLine (userParams, argc - 1, argv + 1);
        std::cout << "=====================================================\n";
        std::cout << userParams;
        std::cout << "=====================================================\n";
    }
    catch (std::runtime_error& error)
    {
        std::cout << error.what() << std::endl;
        return;
    }

    configureFromParameters();
    printStartupMessage();

    while (shouldContinue())
    {
        double dt = getTimestep();

        if (shouldWriteOutput())
        {
            std::string filename = makeOutputFilename();

            PathHelpers::ensureParentDirectoryExists (filename);
            writeOutput (filename);

            std::cout << "n=" << std::setfill ('0') << std::setw (6) << status.simulationIter << " ";
            std::cout << "t=" << std::setw (4) << std::fixed << status.simulationTime << " ";
            std::cout << "dt=" << std::setw (4) << std::scientific << dt << " ";
            std::cout << "output:" << filename << std::endl;

            ++status.outputsWrittenSoFar;
        }

        advance (dt);
        status.simulationTime += dt;
        status.simulationIter += 1;
    }
}

Variant SimulationDriver::getParameter (std::string parameterName) const
{
    return userParams.at (parameterName);
}

SimulationDriver::Status SimulationDriver::getStatus() const
{
    return status;
}

void SimulationDriver::writeAsciiTable (std::vector<std::vector<double>> columns, std::ostream& stream) const
{
    for (int n = 0; n < columns[0].size(); ++n)
    {
        for (int i = 0; i < columns.size(); ++i)
        {
            stream
            << std::scientific
            << std::showpos
            << std::setprecision (10)
            << columns[i][n] << " ";
        }
        stream << std::endl;
    }
}
