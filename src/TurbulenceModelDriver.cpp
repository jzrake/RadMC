#include <fstream>
#include <iomanip>
#include <sstream>
#include "TurbulenceModelDriver.hpp"




// ============================================================================
void TurbulenceModelDriver::makeUserParameters (Variant::NamedValues& params)
{
    params["outdir"] = ".";
    params["tmax"] = 1.0;
    params["kmax"] = 1e4;
    params["bins"] = 128;
    params["cfl"] = 0.5;
    params["cpi"] = 0.1;
    params["urad"] = 1e-2;
    params["lstar"] = 1e-3;
}

void TurbulenceModelDriver::configureFromParameters()
{
    cascade = RichardsonCascade (getParameter ("kmax"), getParameter ("bins"));
    cascade.photonMeanFreePath = getParameter ("lstar");
    cascade.radiativeEnergyDensity = getParameter ("urad");
    cascade.cascadePower = 1.0;
}

void TurbulenceModelDriver::printStartupMessage() const
{
    std::cout << "Photon mean free path: " << cascade.getPhotonMeanFreePathScale() << "\n";
    std::cout << "Viscous scale: " << cascade.getFiducialViscousScale() << "\n";
    std::cout << "Compton power: " << cascade.getFiducialComptonPower() << "\n";
}

double TurbulenceModelDriver::getTimestep() const
{
    double cfl = getParameter ("cfl");
    return cfl * cascade.getShortestTimeScale();
}

void TurbulenceModelDriver::advance (double dt)
{
    cascade.advance (dt);
}

bool TurbulenceModelDriver::shouldContinue() const
{
    SimulationDriver::Status S = getStatus();
    return S.simulationTime < double (getParameter ("tmax"));
}

bool TurbulenceModelDriver::shouldWriteOutput() const
{
    SimulationDriver::Status S = getStatus();
    double timeBetweenOutputs = getParameter ("cpi");
    return S.simulationTime >= timeBetweenOutputs * S.outputsWrittenSoFar - 1e-12;
}

void TurbulenceModelDriver::writeOutput (std::string filename) const
{
    std::vector<std::vector<double>> columns;
    columns.push_back (cascade.powerSpectrum.getDataX());
    columns.push_back (cascade.powerSpectrum.getDataY());
    columns.push_back (cascade.getEddyTurnoverTime());
    columns.push_back (cascade.getDampingTime());

    std::ofstream stream (filename);
    writeAsciiTable (columns, stream);
}

std::string TurbulenceModelDriver::makeOutputFilename() const
{
    SimulationDriver::Status S = getStatus();
    std::ostringstream filenameStream;
    filenameStream << getParameter ("outdir") << "/spectrum.";
    filenameStream << std::setfill ('0') << std::setw (6) << S.outputsWrittenSoFar;
    filenameStream << ".dat";
    return filenameStream.str();;
}
