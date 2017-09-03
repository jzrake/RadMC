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
    using Photon = TurbulentComptonizationModel::Photon;

    m.doc() = "Python module for radiative Monte Carlo calculations";

    pybind11::bind_vector<std::vector<double>> (m, "std::vector<double>", pybind11::module_local (true));
    pybind11::bind_vector<std::vector<FourVector>> (m, "std::vector<FourVector>", pybind11::module_local (true));
    pybind11::bind_vector<std::vector<Photon>> (m, "std::vector<TurbulentComptonizationModel.Photon>", pybind11::module_local (true));

    pybind11::class_<FourVector> fv (m, "FourVector");
    fv.def (pybind11::init<double, double, double, double>());
    fv.def_property_readonly ("t", [] (const FourVector& self) { return self[0]; });
    fv.def_property_readonly ("x", [] (const FourVector& self) { return self[1]; });
    fv.def_property_readonly ("y", [] (const FourVector& self) { return self[2]; });
    fv.def_property_readonly ("z", [] (const FourVector& self) { return self[3]; });

    pybind11::class_<Photon> ph (m, "TurbulentComptonizationModel.Photon");
    ph.def_readwrite ("position", &Photon::position);
    ph.def_readwrite ("momentum", &Photon::momentum);
    ph.def_readwrite ("fluidParcelFourVelocity", &Photon::fluidParcelFourVelocity);

    pybind11::class_<TurbulentComptonizationModel::Config> tcc (m, "TurbulentComptonizationModel.Config");
    tcc.def_readwrite ("outdir", &TurbulentComptonizationModel::Config::outdir);
    tcc.def_readwrite ("tmax", &TurbulentComptonizationModel::Config::tmax);
    tcc.def_readwrite ("cpi", &TurbulentComptonizationModel::Config::cpi);
    tcc.def_readwrite ("tsi", &TurbulentComptonizationModel::Config::tsi);
    tcc.def_readwrite ("theta", &TurbulentComptonizationModel::Config::theta);
    tcc.def_readwrite ("ell_star", &TurbulentComptonizationModel::Config::ell_star);
    tcc.def_readwrite ("nphot", &TurbulentComptonizationModel::Config::nphot);
    tcc.def_readwrite ("nphot_per_mass", &TurbulentComptonizationModel::Config::nphot_per_mass);
    tcc.def_readwrite ("nelec", &TurbulentComptonizationModel::Config::nelec);
    tcc.def_readwrite ("ephot", &TurbulentComptonizationModel::Config::ephot);
    tcc.def_readwrite ("beta_turb", &TurbulentComptonizationModel::Config::beta_turb);

    pybind11::class_<TurbulentComptonizationModel> tcm (m, "TurbulentComptonizationModel");
    tcm.def (pybind11::init<TurbulentComptonizationModel::Config>());
    tcm.def_static ("Config", [] () { return TurbulentComptonizationModel::Config(); });
    tcm.def_static ("Photon", [] () { return TurbulentComptonizationModel::Photon(); });
    tcm.def ("get_cascade_wavenumbers", &TurbulentComptonizationModel::getCascadeWaveNumbers);
    tcm.def ("get_cascade_power_spectrum", &TurbulentComptonizationModel::getCascadePowerSpectrum);
    tcm.def ("get_photons", &TurbulentComptonizationModel::getPhotons);
    tcm.def ("get_photon_energy_spectrum", &TurbulentComptonizationModel::getPhotons);
    tcm.def ("get_electron_temperature", &TurbulentComptonizationModel::getPhotons);
    tcm.def ("get_specific_internal_energy", &TurbulentComptonizationModel::getPhotons);
    tcm.def ("get_specific_kinetic_energy", &TurbulentComptonizationModel::getPhotons);
    tcm.def ("get_specific_radiative_energy", &TurbulentComptonizationModel::getPhotons);
    tcm.def ("get_compton_cooling_time", &TurbulentComptonizationModel::getPhotons);
}

#endif // RADMC_PYTHON
