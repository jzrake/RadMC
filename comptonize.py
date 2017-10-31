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
        energies = np.linspace(min(E), max(E), 1024)

        plt.hist(E, bins=128, histtype='step', normed=True)
        plt.plot(energies, [self.wein_photon_spectrum(temperature, e) for e in energies])
        plt.xlabel(r"$h \nu$")
        plt.ylabel(r"$dN / d\nu$")

    def wein_photon_spectrum(self, kT, nu):
        return 1. / (2 * kT**3) * nu**2 * math.exp(-nu / kT)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    comptonization = ComptonizationTest()
    comptonization.plot_energy_spectrum(1., 200000)
    plt.show()
