import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

from scipy import signal



def smooth(sig):
	win = signal.hann(75)
	return signal.convolve(sig, win, mode='same') / sum(win)



def plot_spectrum(ax, datadir, which):
	fname = os.path.join(datadir, "spectrum.{:06d}.dat".format(which))
	data = np.loadtxt(fname)
	E = data[:,1]
	n = data[:,2] # dN/dE
	#n = smooth(n)
	ax.plot(E, n * E, lw=0.5)

	ax.set_xlabel(r"Photon energy ($m_e c^2$)")
	ax.set_ylabel(r"$dN/d\log\varepsilon$")
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(1e-2, 1e2)
	ax.set_ylim(1e-3, 1e0)



def plot_cascade(ax, datadir, which):
	fname = os.path.join(datadir, "cascade.{:06d}.dat".format(which))
	data = np.loadtxt(fname)
	k = data[:,0] # k
	P = data[:,1] # dP/dk
	ax.plot(k, P, lw=0.5)
	#ax.plot(k, k**(-5./3))

	ax.set_xlabel(r"Inverse scale $\ell^{-1}$")
	ax.set_ylabel(r"$dP/dk$")
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(0.5, 1e3)
	ax.set_ylim(1e-5, 2)



def everything(fig, datadir, which_spec):
	ax_temp = fig.add_subplot(2, 2, 1)
	ax_ener = fig.add_subplot(2, 2, 3)
	ax_spec = fig.add_subplot(2, 2, 2)
	ax_turb = fig.add_subplot(2, 2, 4)

	cpi = 0.02
	fname = os.path.join(datadir, "time-series.dat")

	(simulationTime,
		photonTemperature,
		electronTemperature,
		specificPhotonEnergy,
		specificInternalEnergy,
		specificTurbulentEnergy) = np.loadtxt(fname, unpack=True)

	ax_temp.plot(simulationTime, electronTemperature, label=r'$k T_\pm$')
	ax_temp.plot(simulationTime, photonTemperature, label=r'$k T_\gamma}$')
	ax_temp.legend(loc='best')

	ax_ener.plot(simulationTime, specificPhotonEnergy, label=r'$E_{\gamma}$')
	ax_ener.plot(simulationTime, specificInternalEnergy, label=r'$E_{\rm int}$')
	ax_ener.plot(simulationTime, specificTurbulentEnergy, label=r'$E_{\rm turb}$')
	ax_ener.legend(loc='best')
	ax_ener.set_xlabel("Number of collisions")

	ax_temp.axvline(which_spec * cpi, ymin=0.0, ymax=0.2, ls='--', lw=0.5, c='k')
	ax_ener.axvline(which_spec * cpi, ymin=0.0, ymax=1.0, ls='--', lw=0.5, c='k')

	plot_spectrum(ax_spec, datadir, which_spec)
	plot_cascade(ax_turb, datadir, which_spec)



def make_movie():
	fig = plt.figure(figsize=[10, 6])
	frames = []
	for i in range(1, 150):
		everything(fig, "data/test", i)
		pngname = "{0:04d}.png".format(i)
		frames.append(pngname)
		print(pngname)
		fig.savefig(pngname)
		fig.clf()

	os.system("ffmpeg -i %04d.png -f mp4 -vcodec h264 -pix_fmt yuv420p out.mp4")

	for frame in frames:
		os.remove(frame)


def make_figure():
	fig = plt.figure(figsize=[10, 6])
	everything(fig, "data/test", 20)
	plt.show()


#fig = plt.figure(figsize=[10, 6])
#plot_cascade (fig.add_subplot(1, 1, 1), "data/test", 10)
#plt.show()
#make_figure()
#make_movie()
