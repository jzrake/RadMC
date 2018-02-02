import numpy as np
import matplotlib.pyplot as plt
import radmc



class DissipativeWindPlot(object):
	def __init__(self, fig):
		self.axu = fig.add_subplot(3, 1, 1)
		self.axn = fig.add_subplot(3, 1, 2)
		self.axT = fig.add_subplot(3, 1, 3)

		self.axu.set_xscale('log')
		self.axn.set_xscale('log')
		self.axT.set_xscale('log')
		self.axu.set_yscale('log')
		self.axn.set_yscale('log')
		self.axT.set_yscale('log')

		self.axT.set_xlabel(r'$r / r_0$')
		self.axu.set_ylabel(r'$u = \gamma \beta$')
		self.axn.set_ylabel(r'$\mu_f$')
		self.axT.set_ylabel(r'$\theta = k T / m_e c^2$')

		plt.setp(self.axu.get_xticklabels(), visible=False)
		plt.setp(self.axn.get_xticklabels(), visible=False)

		fig.suptitle("Turbulent Wind Solutions, $n_\gamma / n_p = 10^4$")

	def configure(self, state):
		return state.set_photons_per_baryon(1e4)

	def plot_wind(self, free_enthalpy, heating_rate, label):
		wind = radmc.RelativisticWind()
		wind.set_specific_wind_power(100)
		wind.set_initial_four_velocity(1.0)
		wind.set_initial_free_enthalpy(free_enthalpy)
		wind.set_heating_rate(heating_rate)
		solution = [self.configure(s) for s in wind.integrate_range(1e10)]

		r = [s.r for s in solution]
		u = [s.u for s in solution]
		n = [s.n for s in solution]
		T = [s.temperature() for s in solution]

		self.axu.plot(r, u, label=label)
		self.axT.plot(r, T, label=label)
		if heating_rate > 0:
			self.axn.plot(r, n, label=label)

fig = plt.figure(figsize=[6, 10])
plot = DissipativeWindPlot(fig)
plot.plot_wind(0.00, 0.0, 'ideal')
plot.plot_wind(1e+1, 0.1, r'$\zeta=0.1$')
plot.plot_wind(1e+1, 0.2, r'$\zeta=0.2$')
plot.plot_wind(1e+1, 0.4, r'$\zeta=0.4$')
plot.plot_wind(1e+1, 0.8, r'$\zeta=0.8$')
plot.axu.legend(loc='best')

plt.show()
