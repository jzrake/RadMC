#ifndef SimulationDriver_hpp
#define SimulationDriver_hpp

#include <vector>
#include "Variant.hpp"




class SimulationDriver
{
public:
    struct Status
    {
        Status();
        double simulationTime;
        int simulationIter;
        int outputsWrittenSoFar;
        int timeSeriesSamplesSoFar;
    };

    /**
    Run the simulation, using command line arguments for configuration.
    */
    void run (int argc, const char *argv[]);

    /**
    Get a parameter by name. This will be one of the named values provided by
    makeUserParameters. Raises an exception if the parameter is not found.
    */
    Variant getParameter (std::string parameterName) const;

    /**
    Return a copy of the current simulation status structure.
    */
    Status getStatus() const;

    /**
    Add a named data series that will be tracked over the course of the
    simulation.
    */
    void addTimeSeries (std::string seriesName);

    /**
    Write an ASCII formatted table of the current time series data.
    */
    void writeTimeSeriesData (std::ostream& stream) const;

    /**
    Write an ASCII formatted table of floating point data to the given stream.
    Each column of data is assumed to have the same size.
    */
    void writeAsciiTable (std::vector<std::vector<double>> columns, std::ostream& stream) const;

    virtual void makeUserParameters (Variant::NamedValues&) = 0;
    virtual void configureFromParameters() = 0;
    virtual void printStartupMessage() const = 0;
    virtual double getTimestep() const = 0;
    virtual bool shouldContinue() const = 0;
    virtual void advance (double dt) = 0;
    virtual bool shouldRecordIterationInTimeSeries() const { return false; }
    virtual double getRecordForTimeSeries (std::string) const { return 0.0; }
    virtual bool shouldWriteOutput() const = 0;
    virtual void writeOutput (std::string filename) const = 0;
    virtual std::string makeOutputFilename() const = 0;

private:
    Variant::NamedValues userParams;
    Status status;
    std::vector<std::string> timeSeriesNames;
    std::vector<std::vector<double>> timeSeriesData;
};

#endif
