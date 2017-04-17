#include <iostream>
#include <iomanip>
#include <sstream>
#include "SimulationDriver.hpp"




// ============================================================================
SimulationDriver::Status::Status()
{
    outputsWrittenSoFar = 0;
    timeSeriesSamplesSoFar = 0;
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

        if (shouldRecordIterationInTimeSeries())
        {
            for (int n = 0; n < timeSeriesNames.size(); ++n)
            {
                timeSeriesData[n].push_back (getRecordForTimeSeries (timeSeriesNames[n]));
            }
            ++status.timeSeriesSamplesSoFar;
        }

        if (shouldWriteOutput())
        {
            writeOutput();

            std::cout << "[" << std::setfill ('0') << std::setw (6) << status.simulationIter << "] ";
            std::cout << "t=" << std::setprecision (2) << std::fixed << status.simulationTime << " ";
            std::cout << "dt=" << std::setprecision (2) << std::scientific << dt << "\n";

            ++status.outputsWrittenSoFar;
        }

        advance (dt);
        status.simulationTime += dt;
        status.simulationIter += 1;
    }
}

Variant SimulationDriver::getParameter (std::string parameterName) const
{
    try
    {
       return userParams.at (parameterName);
    }
    catch (std::out_of_range& error)
    {
        std::cerr
        << "[SimulationDriver::getParameter] unknown parameter '"
        << parameterName
        << "'"
        << std::endl;
        throw;
    }
}

SimulationDriver::Status SimulationDriver::getStatus() const
{
    return status;
}

void SimulationDriver::addTimeSeries (std::string seriesName)
{
    timeSeriesNames.push_back (seriesName);
    timeSeriesData.push_back (std::vector<double>());
}

void SimulationDriver::writeTimeSeriesData (std::ostream& stream) const
{
    stream << "#";

    for (auto& seriesName : timeSeriesNames)
    {
        stream << " " << seriesName;
    }
    stream << "\n";

    writeAsciiTable (timeSeriesData, stream);
}

std::string SimulationDriver::makeFilename (std::string directory, std::string base, std::string extension, int number) const
{
    std::ostringstream filenameStream;
    filenameStream << directory << "/" << base;

    if (number >= 0)
    {
        filenameStream << "." << std::setfill ('0') << std::setw (6) << number;
    }
    filenameStream << extension;
    return filenameStream.str();
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
