import argparse
import datetime
import numpy as np
import matplotlib.pyplot as plt
import h5py
import radmc



def make_config():
    config = radmc.StructuredJetModel.Config()
    config.table_resolution_radius = 256
    config.table_resolution_theta = 256
    config.outermost_radius = 1e6
    config.jet_opening_angle = 0.3
    config.jet_polar_boundary = 0.8
    config.jet_structure_exponent = 1.0
    config.specific_wind_power = 50
    config.luminosity_per_steradian = 1e51
    config.inner_radius_cm = 1e8
    config.leptons_per_baryon = 1.0
    config.photons_per_baryon = 1e4
    return config



def report_model(model, config):
    physics = radmc.PhysicsConstants()
    rphot = model.approximate_photosphere(0.0) / config.inner_radius_cm
    mpc2 = physics.gram_to_erg(physics.mp)

    print("jet launch radius .............. {:4e} cm".format(config.inner_radius_cm))
    print("jet photosphere radius ......... {:4e} cm".format(model.approximate_photosphere(0.0)))
    print("baryon number rate ............. {:4e} 1/s/Sr".format(config.luminosity_per_steradian / config.specific_wind_power / mpc2))
    print("on-axis luminosity per Sr ...... {:4e} erg/s/Sr".format(config.luminosity_per_steradian))
    print("total luminosity ............... {:4e} erg/s".format(model.total_luminosity()))
    print("lag time on-axis ............... {:4e} s".format(model.approximate_lag_time(0.0)))
    print("lag time @ theta=0.5 ........... {:4e} s".format(model.approximate_lag_time(0.5)))
    print("fluid travel to photosphere .... {:4e} s".format(model.fluid_propagation_time_to_radius(rphot, 0.0)))
    print("light travel to photosphere .... {:4e} s".format(rphot * config.inner_radius_cm / physics.c))



class StructuredJetModel(object):
    def __init__(self, filename, restart=False):

        if restart:
            self.file = h5py.File(filename, 'a')
            self.photon_number = len(self.file['time'])
            config = self.read_config_from_file()

        else:
            config = make_config()
            self.file = h5py.File(filename, 'w')
            self.photon_number = 0
            self.write_config_to_file(config)

            for key in ['time', 'radius', 'theta', 'energy', 'lag']:
                if key not in self.file:
                    self.file.create_dataset(key, (0,), maxshape=(None,))

        self.config = config
        self.model = radmc.StructuredJetModel(config)
        self.photon_cache = []
        report_model(self.model, config)

    def run_photon(self):
        q = self.model.sample_theta(0.5)
        r = self.model.approximate_photosphere(q) / self.config.inner_radius_cm * 0.01
        p = self.model.generate_photon(r, q)
        n = 0

        while (p.position.radius < self.config.outermost_radius * 0.1
            and p.position.theta < self.config.jet_polar_boundary):
            p = self.model.step_photon(p)
            n += 1

        self.photon_cache.append(p)
        self.photon_number += 1
        r0 = self.config.inner_radius_cm

        print("photon {} @ q : {:.2f} -> {:.2f} Nscat = {}, final r and E are {:.2e}, {:.2e}, lag = {:.2e}s"
            .format(self.photon_number, q, p.position.theta, n, p.position.radius, p.momentum.t, p.lag_time(r0)))

    def purge_photons(self):
        n0 = len(self.file['time'])
        n1 = n0 + len(self.photon_cache)
        r0 = self.config.inner_radius_cm

        entries = dict()
        entries['time'  ] = [p.position.t for p in self.photon_cache]
        entries['radius'] = [p.position.radius for p in self.photon_cache]
        entries['theta' ] = [p.momentum.theta for p in self.photon_cache]
        entries['energy'] = [p.momentum.t for p in self.photon_cache]
        entries['lag'   ] = [p.lag_time(r0) for p in self.photon_cache]

        for key in entries:
            self.file[key].resize((n1,))
            self.file[key][n0:] = entries[key]

        self.photon_cache = []
        print("purged photons, now {} in file".format(n1))

    def write_config_to_file(self, config):
        config_group = self.file.create_group('config')
        for key in dir(config):
            if key.startswith('__'): continue
            config_group[key] = getattr(config, key)

    def read_config_from_file(self):
        config_group = self.file['config']
        config = radmc.StructuredJetModel.Config()
        for key in dir(config):
            if key.startswith('__'): continue
            setattr(config, key, config_group[key].value)
        return config



class WindProfile(object):

    def __init__(self):
        config = make_config()
        self.config = config
        self.model = radmc.StructuredJetModel(config)
        report_model(self.model, config)

    def plot_wind_profile(self):
        fig = plt.figure(figsize=[9, 6])
        ax1 = fig.add_subplot(1, 1, 1)

        r = np.logspace(0, 6, 128)
        cm = self.config.inner_radius_cm

        for theta in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
            states = [self.model.sample_wind_spherical(ri, theta) for ri in r]
            rphot = self.model.approximate_photosphere(theta)
            uphot = self.model.sample_wind_spherical(rphot / cm, theta).u
            lag = self.model.approximate_lag_time(theta)
            ax1.plot(r * cm, [s.u for s in states])
            ax1.scatter([rphot], [uphot], marker='*', c='y', s=200)
            ax1.text(rphot * 0.75, uphot * 1.25, r"$\delta_t = {:.1e}s$".format(lag))
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r"radius, cm")
        ax1.set_ylabel(r"four velocity (c)")
        plt.show()



class Analysis(object):

    def __init__(self, filename):
        self.file = h5py.File(filename, 'r')

    def make_log_histogram(self, E):
        dy, x = np.histogram(E, bins=np.logspace(np.log10(min(E)), np.log10(max(E)), 100))
        dx = (x[1:] - x[:-1])
        xc = 0.5 * (x[1:] + x[:-1])
        return xc, dy / dx

    def make_lin_histogram(self, E):
        dy, x = np.histogram(E, bins=np.linspace(min(E), max(E), 100))
        dx = (x[1:] - x[:-1])
        xc = 0.5 * (x[1:] + x[:-1])
        return xc, dy / dx

    def load_photons_for_viewing_angles(self, theta0, theta1):
        T = self.file['time'][...]
        R = self.file['radius'][...]
        Q = self.file['theta'][...]
        E = self.file['energy'][...]
        L = self.file['lag'][...]
        I = np.where((theta0 < Q) * (Q <= theta1))# * (L < 10.0))
        return T[I], R[I], Q[I], E[I], L[I]

    def plot_stats(self):
        fig = plt.figure(figsize=[10, 8])
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)

        angles = [0.0, 0.2, 0.4, 0.6]

        for i in range(len(angles) - 1):
            T, R, Q, E, L = self.load_photons_for_viewing_angles(angles[i], angles[i + 1])

            E0, NE = self.make_log_histogram(E)
            Q0, NQ = self.make_log_histogram(Q)

            # if i == 0: ax1.loglog(E0, 1e6 * E0**1.8)
            ax1.loglog(E0, NE * E0**2)
            ax2.loglog(Q0, NQ / (2 * np.pi * np.sin(Q0)))
            ax3.hist(L, histtype='step', bins=128, log=False)
            ax4.scatter(L, Q, s=.1)

        ax1.set_ylabel(r"$E^2 dN/dE$")
        ax2.set_ylabel(r"$dN/d\Omega$")
        ax3.set_ylabel(r"photon flux")
        ax3.set_xlabel(r"$t - t_c$")
        ax4.set_xlabel(r"lag time")
        ax4.set_ylabel(r"observer angle")
        #ax3.set_xlim(0.0, 10.0)
        #ax3.set_ylim(1e-1, 100.0)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=['run', 'plot', 'plot_wind'])
    parser.add_argument("--photons", type=int, default=1000)
    parser.add_argument("--restart", action='store_true')
    parser.add_argument("--file", default="structured_jet.h5")
    args = parser.parse_args()

    if args.command == 'run':
        radmc.seed(str(datetime.datetime.now()))
        model = StructuredJetModel(args.file, restart=args.restart)

        while model.photon_number < args.photons:

            try:
                model.run_photon()
            except RuntimeError as e:
                print("Bad photon:", e)

            if len(model.photon_cache) >= 10:
                model.purge_photons()

    elif args.command == 'plot':
        analysis = Analysis(args.file)
        analysis.plot_stats()
        plt.show()

    elif args.command == 'plot_wind':
        WindProfile().plot_wind_profile()






# def make_estimates(theta):
#     theta_jet = 0.1
#     viewing_angle = 0.24

#     physics = radmc.PhysicsConstants()
#     Etot = 2e51 # erg
#     L0 = 5e52 # erg / s / Sr
#     eta0 = 1000
#     f0 = L0 / eta0 # erg / s / Sr

#     f = f0 * np.exp(-(theta / theta_jet)**2)
#     eta = eta0 * np.exp(-(theta / theta_jet)**2)

#     Ndot = f / (physics.mp * physics.c * physics.c)
#     rphot = Ndot * physics.st / (2 * eta**2 * physics.c)
#     tlag = rphot / (2 * eta**2 * physics.c)

#     print("Ndot .......... {:.2e} 1/s/Sr".format(Ndot))
#     print("rphot ......... {:.2e} cm".format(rphot))
#     print("tlag .......... {:.2e} s".format(tlag))
#     print("L ............. {:.2e} s".format(eta * f))

# make_estimates(0.24)
# exit()

