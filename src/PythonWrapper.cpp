#ifdef RADMC_PYTHON
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/numpy.h"

#include "PythonWrapper.hpp"
#include "FourVector.hpp"
#include "RichardsonCascade.hpp"
#include "TurbulentComptonizationModel.hpp"





PYBIND11_MODULE (radmc, m)
{
    using TCM = TurbulentComptonizationModel;
    namespace py = pybind11;

    m.doc() = "Python module for radiative Monte Carlo calculations";

    py::bind_vector<std::vector<double>> (m, "std::vector<double>", py::module_local (true));
    py::bind_vector<std::vector<FourVector>> (m, "std::vector<FourVector>", py::module_local (true));
    py::bind_vector<std::vector<TCM::Photon>> (m, "std::vector<TurbulentComptonizationModel.Photon>", py::module_local (true));

    py::class_<FourVector> fv (m, "FourVector");
    fv.def (py::init<double, double, double, double>());
    fv.def_property_readonly ("t", [] (const FourVector& self) { return self[0]; });
    fv.def_property_readonly ("x", [] (const FourVector& self) { return self[1]; });
    fv.def_property_readonly ("y", [] (const FourVector& self) { return self[2]; });
    fv.def_property_readonly ("z", [] (const FourVector& self) { return self[3]; });

    auto tcm = py::class_<TCM> (m, "TurbulentComptonizationModel", py::dynamic_attr())
    .def (py::init<TCM::Config>())
    .def ("advance", &TCM::advance)
    .def ("get_timestep", &TCM::getTimestep)
    .def ("get_photons", &TCM::getPhotons)
    .def ("get_cascade_wavenumber_bins", &TCM::getCascadeWaveNumberBins)
    .def ("get_cascade_power_spectrum", &TCM::getCascadePowerSpectrum)
    .def ("get_photon_energy_bins", &TCM::getPhotonEnergyBins)
    .def ("get_photon_spectrum", &TCM::getPhotonSpectrum)
    .def ("get_electron_temperature", &TCM::getElectronTemperature)
    .def ("get_photon_temperature", &TCM::getPhotonTemperature)
    .def ("get_effective_wave_temperature", &TCM::getEffectiveWaveTemperature)
    .def ("get_compton_cooling_time", &TCM::getComptonCoolingTime)
    .def ("get_specific_internal_energy", &TCM::getSpecificInternalEnergy)
    .def ("get_specific_kinetic_energy", &TCM::getSpecificKineticEnergy)
    .def ("get_specific_photon_energy", &TCM::getSpecificPhotonEnergy)
    .def ("get_eddy_velocity_at_scale", &TCM::getEddyVelocityAtScale);

    py::class_<TCM::Photon> (tcm, "Photon")
    .def (py::init<>())
    .def_readwrite ("position", &TCM::Photon::position)
    .def_readwrite ("momentum", &TCM::Photon::momentum)
    .def_readwrite ("fluidParcelFourVelocity", &TCM::Photon::fluidParcelFourVelocity);

    py::class_<TCM::Config> (tcm, "Config")
    .def (py::init<>())
    .def_readwrite ("theta", &TCM::Config::theta)
    .def_readwrite ("ell_star", &TCM::Config::ell_star)
    .def_readwrite ("nphot", &TCM::Config::nphot)
    .def_readwrite ("nphot_per_mass", &TCM::Config::nphot_per_mass)
    .def_readwrite ("nelec", &TCM::Config::nelec)
    .def_readwrite ("ephot", &TCM::Config::ephot)
    .def_readwrite ("beta_turb", &TCM::Config::beta_turb);

    py::class_<TCM::IterationReport> (tcm, "IterationReport")
    .def (py::init<>())
    .def_readonly ("time_after_iteration", &TCM::IterationReport::timeAfterIteration)
    .def_readonly ("mean_scatterings_per_photon", &TCM::IterationReport::meanScatteringsPerPhoton)
    .def_readonly ("mean_scattering_angle_with_bulk", &TCM::IterationReport::meanScatteringAngleWithBulk)
    .def_readonly ("mean_scattering_angle_in_parcel", &TCM::IterationReport::meanScatteringAngleInParcel)
    .def_readonly ("viscous_power_based_on_photons", &TCM::IterationReport::viscousPowerBasedOnPhotons)
    .def_readonly ("viscous_power_based_on_cascade", &TCM::IterationReport::viscousPowerBasedOnCascade);
}

#endif // RADMC_PYTHON
