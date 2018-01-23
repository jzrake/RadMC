#ifdef RADMC_PYTHON
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/numpy.h"

#include "PythonWrapper.hpp"
#include "FourVector.hpp"
#include "RichardsonCascade.hpp"
#include "TurbulentComptonizationModel.hpp"
#include "StructuredJetModel.hpp"
#include "ScatteringOperations.hpp"
#include "PhysicsConstants.hpp"




PYBIND11_MODULE (radmc, m)
{
    namespace py = pybind11;
    m.doc() = "Python module for radiative Monte Carlo calculations";


    // ========================================================================
    // Generic data structures
    // ========================================================================
    py::bind_vector<std::vector<double>> (m, "std::vector<double>", py::module_local (true));
    py::bind_vector<std::vector<FourVector>> (m, "std::vector<FourVector>", py::module_local (true));

    m.def ("seed", RandomVariable::seed);

    py::class_<FourVector> (m, "FourVector")
    .def (py::init<double, double, double, double>())
    .def_property ("t", [] (const FourVector& self) { return self[0]; }, [] (FourVector& self, double a) { self[0] = a; })
    .def_property ("x", [] (const FourVector& self) { return self[1]; }, [] (FourVector& self, double a) { self[1] = a; })
    .def_property ("y", [] (const FourVector& self) { return self[2]; }, [] (FourVector& self, double a) { self[2] = a; })
    .def_property ("z", [] (const FourVector& self) { return self[3]; }, [] (FourVector& self, double a) { self[3] = a; })
    .def_property_readonly ("radius", [] (const FourVector& self) { return self.radius(); })
    .def_property_readonly ("theta", [] (const FourVector& self) { return self.theta(); });

    py::class_<PhysicsConstants> (m, "PhysicsConstants")
    .def (py::init<>())
    .def_readonly ("c", &PhysicsConstants::c)
    .def_readonly ("k", &PhysicsConstants::k)
    .def_readonly ("h", &PhysicsConstants::h)
    .def_readonly ("G", &PhysicsConstants::G)
    .def_readonly ("e", &PhysicsConstants::e)
    .def_readonly ("me", &PhysicsConstants::me)
    .def_readonly ("mp", &PhysicsConstants::mp)
    .def_readonly ("st", &PhysicsConstants::st)
    .def_readonly ("pc", &PhysicsConstants::pc)
    .def ("gram_to_erg", &PhysicsConstants::gramToErg)
    .def ("gram_to_MeV", &PhysicsConstants::gramToMeV);

    py::class_<ScatteringOperations> (m, "ScatteringOperations")
    .def (py::init<>())
    .def ("sample_scattered_particles", &ScatteringOperations::sampleScatteredParticles)
    .def ("sample_scattered_particles_in_frame", &ScatteringOperations::sampleScatteredParticlesInFrame)
    .def ("compton_scatter", &ScatteringOperations::comptonScatter)
    .def ("comptonize", &ScatteringOperations::comptonize);


    // ========================================================================
    // TurbulentComptonizationModel
    // ========================================================================
    using TCM = TurbulentComptonizationModel;
    py::bind_vector<std::vector<TCM::Photon>> (m, "std::vector<TurbulentComptonizationModel.Photon>", py::module_local (true));

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
    .def ("get_eddy_velocity_at_scale", &TCM::getEddyVelocityAtScale)
    .def ("get_average_compton_y", &TCM::getAverageComptonY);

    py::enum_<TCM::ElectronTemperatureMode> (tcm, "ElectronTemperatureMode")
    .value("Consistent", TCM::ElectronTemperatureMode::Consistent)
    .value("Cold", TCM::ElectronTemperatureMode::Cold)
    .value("LockPhoton", TCM::ElectronTemperatureMode::LockPhoton)
    .value("LockUser", TCM::ElectronTemperatureMode::LockUser);

    py::class_<TCM::Photon> (tcm, "Photon")
    .def (py::init<>())
    .def_readwrite ("position", &TCM::Photon::position)
    .def_readwrite ("momentum", &TCM::Photon::momentum)
    .def_readwrite ("fluidParcelFourVelocity", &TCM::Photon::fluidParcelFourVelocity);

    py::class_<TCM::Config> (tcm, "Config")
    .def (py::init<>())
    .def_readwrite ("electron_temperature_mode", &TCM::Config::electron_temperature_mode)
    .def_readwrite ("disable_cascade_model", &TCM::Config::disable_cascade_model)
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


    // ========================================================================
    // RelativisticWind
    // ========================================================================
    using SRW = RelativisticWind;
    auto srw = py::class_<SRW> (m, "RelativisticWind", py::dynamic_attr());
    srw.def (py::init<>());
    srw.def ("integrate", &SRW::integrate);
    srw.def ("set_specific_wind_power", &SRW::setSpecificWindPower);
    srw.def ("set_initial_four_velocity", &SRW::setInitialFourVelocity);
    srw.def ("set_entropy_production_rate", &SRW::setEntropyProductionRate);

    py::class_<SRW::WindState> (srw, "WindState")
    .def ("set_luminosity_per_steradian", &SRW::WindState::setLuminosityPerSteradian)
    .def ("set_inner_radius_cm", &SRW::WindState::setInnerRadiusCm)
    .def ("set_leptons_per_baryon", &SRW::WindState::setLeptonsPerBaryon)
    .def ("set_photons_per_baryon", &SRW::WindState::setPhotonsPerBaryon)
    .def ("set_propagation_angle", &SRW::WindState::setPropagationAngle)
    .def ("temperature", &SRW::WindState::temperature)
    .def ("proper_number_density", &SRW::WindState::properNumberDensity)
    .def ("blackbody_photons_per_proton", &SRW::WindState::blackbodyPhotonsPerProton)
    .def ("thomson_mean_free_path", &SRW::WindState::thomsonMeanFreePath)
    .def ("four_velocity", &SRW::WindState::fourVelocity)
    .def_readonly ("r", &SRW::WindState::r)
    .def_readonly ("u", &SRW::WindState::u)
    .def_readonly ("g", &SRW::WindState::g)
    .def_readonly ("m", &SRW::WindState::m)
    .def_readonly ("p", &SRW::WindState::p)
    .def_readonly ("d", &SRW::WindState::d)
    .def_readonly ("s", &SRW::WindState::s);

    // ========================================================================
    // StructuredJetModel
    // ========================================================================
    using SJM = StructuredJetModel;
    py::bind_vector<std::vector<SJM::Photon>> (m, "std::vector<StructuredJetModel.Photon>", py::module_local (true));

    auto sjm = py::class_<SJM> (m, "StructuredJetModel", py::dynamic_attr())
    .def (py::init<SJM::Config>())
    .def ("sample_theta", &SJM::sampleTheta)
    .def ("approximate_photosphere", &SJM::approximatePhotosphere)
    .def ("approximate_lag_time", &SJM::approximateLagTime)
    .def ("angular_luminosity", &SJM::angularLuminosity)
    .def ("total_luminosity", &SJM::totalLuminosity)
    .def ("fluid_propagation_time_to_radius", &SJM::fluidPropagationTimeToRadius)
    .def ("generate_photon", &SJM::generatePhoton)
    .def ("step_photon", &SJM::stepPhoton)
    .def ("sample_wind", &SJM::sampleWind)
    .def ("sample_wind_spherical", &SJM::sampleWindSpherical);

    py::class_<SJM::Config> (sjm, "Config")
    .def (py::init<>())
    .def_readwrite ("table_resolution_radius", &SJM::Config::tableResolutionRadius)
    .def_readwrite ("table_resolution_theta", &SJM::Config::tableResolutionTheta)
    .def_readwrite ("outermost_radius", &SJM::Config::outermostRadius)
    .def_readwrite ("jet_opening_angle", &SJM::Config::jetOpeningAngle)
    .def_readwrite ("jet_polar_boundary", &SJM::Config::jetPolarBoundary)
    .def_readwrite ("jet_structure_exponent", &SJM::Config::jetStructureExponent)
    .def_readwrite ("specific_wind_power", &SJM::Config::specificWindPower)
    .def_readwrite ("luminosity_per_steradian", &SJM::Config::luminosityPerSteradian)
    .def_readwrite ("inner_radius_cm", &SJM::Config::innerRadiusCm)
    .def_readwrite ("leptons_per_baryon", &SJM::Config::leptonsPerBaryon)
    .def_readwrite ("photons_per_baryon", &SJM::Config::photonsPerBaryon);

    py::class_<SJM::Photon> (sjm, "Photon")
    .def (py::init<>())
    .def_readwrite ("position", &SJM::Photon::position)
    .def_readwrite ("momentum", &SJM::Photon::momentum)
    .def ("lag_time", &SJM::Photon::lagTime);
}

#endif // RADMC_PYTHON
