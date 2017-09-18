import numpy as np
import h5py
import radmc
import TaskScheduler




class PrintIterationMessage(TaskScheduler.Task):
	def __init__(self, model):
		self.model = model

	def run(self, status, repetition):
		report = self.model.report

		print("[{0:08d}]: t={1:.2e} sc/ph={2:.1e}".format(
			status.simulation_iter,
			status.simulation_time,
			report.mean_scatterings_per_photon))

	def get_recurrence(self):
	    return TaskScheduler.Recurrence(0.0, 0.0, 100)




class CollectTimeSeriesData(TaskScheduler.Task):
	def __init__(self, model):
		self.model = model
		self.h5_time_series = model.h5_file.create_group('time_series')
		for key in self.sample():
			self.h5_time_series.create_dataset(key, (0,), maxshape=(None,))

	def run(self, status, repetition):
		for key, value in self.sample(status.simulation_time).iteritems():
			self.h5_time_series[key].resize((repetition + 1,))
			self.h5_time_series[key][-1] = value

	def get_recurrence(self):
		return TaskScheduler.Recurrence(0.01)

	def sample(self, time=0.0):
		return dict(
			time = time,
			phot_kT                    = self.model.get_photon_temperature(),
			elec_kT                    = self.model.get_electron_temperature(),
			wave_kT                    = self.model.get_effective_wave_temperature(),
			compton_cooling_time       = self.model.get_compton_cooling_time(),
			specific_kinetic_energy    = self.model.get_specific_kinetic_energy(),
			specific_internal_energy   = self.model.get_specific_internal_energy(),
			specific_photon_energy     = self.model.get_specific_photon_energy())




class CollectSpectralData(TaskScheduler.Task):
	def __init__(self, model):
		self.model = model
		self.h5_spectra = model.h5_file.create_group('spectra')

	def run(self, status, repetition):
		itr = self.h5_spectra.create_group('{0:06d}'.format(repetition))
		itr['photon_E']  = self.model.get_photon_energy_bins()
		itr['photon_N']  = self.model.get_photon_spectrum()
		itr['cascade_k'] = self.model.get_cascade_wavenumber_bins()
		itr['cascade_P'] = self.model.get_cascade_power_spectrum()
		itr['time'] = status.simulation_time

	def get_recurrence(self):
		return TaskScheduler.Recurrence(0.2)




cfg = radmc.TurbulentComptonizationModel.Config()
cfg.nphot          = 4
cfg.nphot_per_mass = 1e0
cfg.ephot          = 1e-6
cfg.ell_star       = 0.1
cfg.beta_turb      = 0.1
cfg.theta          = 1e-6


scheduler = TaskScheduler.TaskScheduler()
status = TaskScheduler.Status()
model = radmc.TurbulentComptonizationModel(cfg)
model.report = radmc.TurbulentComptonizationModel.IterationReport()
model.h5_file = h5py.File('radmc.h5', 'w')

h5_cfg = model.h5_file.create_group('config')
for k in dir(cfg):
	if not k.startswith('__'):
		h5_cfg[k] = getattr(cfg, k)


scheduler.schedule(PrintIterationMessage(model))
scheduler.schedule(CollectTimeSeriesData(model))
scheduler.schedule(CollectSpectralData(model))


#viscous_power_based_on_photons = 0.0
#viscous_power_based_on_cascade = 0.0

while status.simulation_time < 80.0:

	scheduler.dispatch(status)
	dt = model.get_timestep()
	model.report = model.advance(dt)
	status.step(dt)

	#viscous_power_based_on_photons += model.report.viscous_power_based_on_photons
	#viscous_power_based_on_cascade += model.report.viscous_power_based_on_cascade
	#print(model.report.viscous_power_based_on_photons, viscous_power_based_on_photons, viscous_power_based_on_cascade)

