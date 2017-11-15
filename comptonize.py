import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import radmc



class ComptonizationTest(object):
    def __init__(self):
        pass

    def plot_energy_spectrum(self, temperature, num_scatterings=10000):
        ops = radmc.ScatteringOperations()
        path = ops.comptonize(temperature, num_scatterings)

        E = [p.t for p in path]

        bins = np.logspace(np.log10(min(E)), np.log10(max(E)), 100)
        ax1 = plt.figure().add_subplot(1, 1, 1)
        ax1.hist(E, bins=bins, histtype='step', normed=True)
        ax1.plot(bins, [self.wien_photon_spectrum(temperature, e) for e in bins])
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r"$h \nu$")
        ax1.set_ylabel(r"$dN / d\nu$")

    def wien_photon_spectrum(self, kT, nu):
        return 1. / (2 * kT**3) * nu**2 * math.exp(-nu / kT)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    comptonization = ComptonizationTest()
    comptonization.plot_energy_spectrum(1., 200000)
    plt.show()
