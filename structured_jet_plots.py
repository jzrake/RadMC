import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt
import radmc


def make_log_histogram(E):
	dy, x = np.histogram(E, bins=np.logspace(np.log10(min(E)), np.log10(max(E)) + 1, 100))
	return 0.5 * (x[1:] + x[:-1]), dy / (x[1:] - x[:-1])


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("file")
	args = parser.parse_args()

	file = h5py.File(args.file, 'r')

	E = []
	T = []
	R = []
	Q = []

	for key in file['photon_paths']:
		E.append(file['photon_paths'][key]['energy'].value)
		#T.append(file['photon_paths'][key]['time'].value)
		#R.append(file['photon_paths'][key]['radius'].value)
		Q.append(file['photon_paths'][key]['theta'].value)

	fig = plt.figure()
	ax1 = fig.add_subplot(2, 1, 1)
	ax2 = fig.add_subplot(2, 1, 2)

	E, NE = make_log_histogram(E)
	ax1.loglog(E, NE * E)
	ax1.loglog(E, 1e3 * E)
	ax1.set_ylabel(r"$E dN/dE$")

	Q, NQ = make_log_histogram(Q)
	ax2.loglog(Q, NQ)
	ax2.set_ylabel(r"$dN/d\theta$")

	plt.show()

