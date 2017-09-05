import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation



class CascadePlot(object):
    def __init__(self, ax, h5_file):
        self.ax = ax
        self.h5_file = h5_file
        self.line, = ax.plot([], [], '--')

    def initialize_axes(self):
        ax = self.ax
        ax.set_xlabel(r"Inverse scale $\ell^{-1}$")
        ax.set_ylabel(r"$dP/dk$")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5, 1e3)
        ax.set_ylim(1e-7, 2)
        return self.line,

    def load_frame(self, frame):
    	grp = self.h5_file['spectra'].values()[frame]
        k = grp['cascade_k']
        P = grp['cascade_P']
        self.line.set_xdata(k)
        self.line.set_ydata(P)
        return self.line,



class SpectrumPlot(object):
    def __init__(self, ax, h5_file):
        self.ax = ax
        self.h5_file = h5_file
        self.line, = ax.plot([], [], '-', lw=0.5, c='k')

    def initialize_axes(self):
        ax = self.ax
        ax.set_xlabel(r"Photon energy ($m_e c^2$)")
        ax.set_ylabel(r"$dN/d\log\varepsilon$")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-6, 1e0)
        ax.set_ylim(1e-3, 1e0)
        return self.line,

    def load_frame(self, frame):
		grp = self.h5_file['spectra'].values()[frame]
		E = grp['photon_E'][...]
		n = grp['photon_N'][...]
		self.line.set_xdata(E)
		self.line.set_ydata(E * n)
		return self.line,



class TimeSeriesTemperatures(object):
    def __init__(self, ax, h5_file):
        self.ax = ax
        self.h5_file = h5_file
        self.simulationTime          = h5_file['time_series']['time'][...]
        self.photonTemperature       = h5_file['time_series']['phot_kT'][...]
        self.electronTemperature     = h5_file['time_series']['elec_kT'][...]

        self.kT_electron,  = self.ax.plot([], [], label=r'$k T_\pm$')
        self.kT_photon,    = self.ax.plot([], [], label=r'$k T_\gamma}$')
        self.vertical_bar, = self.ax.plot([], [], ls='--', lw=0.5, c='k')

    def initialize_axes(self):
        self.kT_electron.set_data(self.simulationTime, self.electronTemperature)
        self.kT_photon.set_data(self.simulationTime, self.photonTemperature)

        self.ax.relim()
        self.ax.autoscale(enable=True)
        self.ax.set_xlabel("Time (eddy turnover)")

        return self.kT_electron, self.kT_photon, self.vertical_bar

    def load_frame(self, frame):
        t = self.h5_file['spectra'].values()[frame]['time'].value

        self.vertical_bar.set_xdata([t, t])
        self.vertical_bar.set_ydata(self.ax.get_ylim())

        return self.kT_electron, self.kT_photon, self.vertical_bar



class TimeSeriesEnergies(object):
    def __init__(self, ax, h5_file):
        self.ax = ax
        self.h5_file = h5_file
        self.simulationTime          = h5_file['time_series']['time'][...]
        self.specificKineticEnergy   = h5_file['time_series']['specific_kinetic_energy'][...]
        self.specificInternalEnergy  = h5_file['time_series']['specific_internal_energy'][...]
        self.specificPhotonEnergy    = h5_file['time_series']['specific_photon_energy'][...]

        self.Eint, = self.ax.plot([], [], label=r'$E_{\rm int}$')
        self.Erad, = self.ax.plot([], [], label=r'$E_\gamma$')
        self.Ekin, = self.ax.plot([], [], label=r'$E_{\rm turb}$')
        self.vertical_bar, = self.ax.plot([], [], ls='--', lw=0.5, c='k')

    def initialize_axes(self):
        self.Eint.set_data(self.simulationTime, self.specificInternalEnergy)
        self.Erad.set_data(self.simulationTime, self.specificPhotonEnergy)
        self.Ekin.set_data(self.simulationTime, self.specificKineticEnergy)

        self.ax.relim()
        self.ax.autoscale(enable=True)
        self.ax.set_xlabel("Time (eddy turnover)")
        self.ax.legend(loc='best')

        return self.Eint, self.Erad, self.Ekin, self.vertical_bar

    def load_frame(self, frame):
        t = self.h5_file['spectra'].values()[frame]['time'].value

        self.vertical_bar.set_xdata([t, t])
        self.vertical_bar.set_ydata(self.ax.get_ylim())

        return self.Eint, self.Erad, self.Ekin, self.vertical_bar



def composite_load_frame(*args):
    def f(frame):
        artists = []
        for arg in args:
            artists += arg.load_frame(frame)
        return artists
    return f



def composite_initialize(*args):
    def f():
        artists = []
        for arg in args:
            artists += arg.initialize_axes()
        return artists
    return f


h5_file = h5py.File('radmc.h5', 'r')
fig = plt.figure(figsize=[10,6])
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)
plot1 = TimeSeriesEnergies(ax1, h5_file)
plot2 = CascadePlot(ax2, h5_file)
plot3 = TimeSeriesTemperatures(ax3, h5_file)
plot4 = SpectrumPlot(ax4, h5_file)

plots = [plot1, plot2, plot3, plot4]

animation = FuncAnimation(
    fig,
    composite_load_frame(*plots),
    init_func=composite_initialize(*plots),
    frames=range(0, len(h5_file['spectra'].keys()), 1),
    interval=10,
    blit=True)

# for plot in plots:
# 	plot.initialize_axes()
# 	plot.load_frame(99)

plt.show()


