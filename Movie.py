import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

cpi = 0.1


class CascadePlot(object):
    def __init__(self, ax):
        self.ax = ax
        self.line, = ax.plot([], [], '--')
        self.datadir = "data/test"

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
        fname = os.path.join(self.datadir, "cascade.{:06d}.dat".format(frame))
        data = np.loadtxt(fname)
        k = data[:,0] # k
        P = data[:,1] # dP/dk
        self.line.set_xdata(k)
        self.line.set_ydata(P)
        return self.line,



class SpectrumPlot(object):
    def __init__(self, ax):
        self.ax = ax
        self.line, = ax.plot([], [], '-', lw=0.5, c='k')
        self.datadir = "data/test"

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
        fname = os.path.join(self.datadir, "spectrum.{:06d}.dat".format(frame))
        data = np.loadtxt(fname)
        E = data[:,1]
        n = data[:,2] # dN/dE
        self.line.set_xdata(E)
        self.line.set_ydata(E * n)
        return self.line,



class TimeSeriesTemperatures(object):
    def __init__(self, ax):
        self.ax = ax
        self.datadir = "data/test"

        fname = os.path.join(self.datadir, "time-series.dat")
        data = np.loadtxt(fname, unpack=True)

        self.simulationTime          = data[0]
        self.photonTemperature       = data[1]
        self.electronTemperature     = data[2]

    def initialize_axes(self):
        self.kT_electron,  = self.ax.plot([], [], label=r'$k T_\pm$')
        self.kT_photon,    = self.ax.plot([], [], label=r'$k T_\gamma}$')
        self.vertical_bar, = self.ax.plot([], [], ls='--', lw=0.5, c='k')

        self.kT_electron.set_data(self.simulationTime, self.electronTemperature)
        self.kT_photon.set_data(self.simulationTime, self.photonTemperature)

        self.ax.relim()
        self.ax.autoscale(enable=True)
        self.ax.set_xlabel("Time (eddy turnover)")

        return self.kT_electron, self.kT_photon, self.vertical_bar

    def load_frame(self, frame):
        t = frame * cpi

        self.vertical_bar.set_xdata([t, t])
        self.vertical_bar.set_ydata(self.ax.get_ylim())

        return self.kT_electron, self.kT_photon, self.vertical_bar



class TimeSeriesEnergies(object):
    def __init__(self, ax):
        self.ax = ax
        self.datadir = "data/test"

        fname = os.path.join(self.datadir, "time-series.dat")
        data = np.loadtxt(fname, unpack=True)

        self.simulationTime          = data[0]
        self.specificPhotonEnergy    = data[3]
        self.specificInternalEnergy  = data[4]
        self.specificTurbulentEnergy = data[5]

    def initialize_axes(self):
        self.Eint, = self.ax.plot([], [], label=r'$E_{\rm int}$')
        self.Erad, = self.ax.plot([], [], label=r'$E_\gamma$')
        self.Ekin, = self.ax.plot([], [], label=r'$E_{\rm turb}$')
        self.vertical_bar, = self.ax.plot([], [], ls='--', lw=0.5, c='k')

        self.Eint.set_data(self.simulationTime, self.specificInternalEnergy)
        self.Erad.set_data(self.simulationTime, self.specificPhotonEnergy)
        self.Ekin.set_data(self.simulationTime, self.specificTurbulentEnergy)

        self.ax.relim()
        self.ax.autoscale(enable=True)
        self.ax.set_xlabel("Time (eddy turnover)")

        return self.Eint, self.Erad, self.Ekin, self.vertical_bar

    def load_frame(self, frame):
        t = frame * cpi

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



fig = plt.figure(figsize=[10,6])
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)
plot1 = TimeSeriesEnergies(ax1)
plot2 = CascadePlot(ax2)
plot3 = TimeSeriesTemperatures(ax3)
plot4 = SpectrumPlot(ax4)

animation = FuncAnimation(
    fig,
    composite_load_frame(plot1, plot2, plot3, plot4),
    init_func=composite_initialize(plot1, plot2, plot3, plot4),
    frames=range(0, 770, 3),
    interval=10,
    blit=True)

plt.show()


