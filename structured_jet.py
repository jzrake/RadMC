#!/usr/bin/env python3

import argparse
import datetime
import numpy as np
import matplotlib.pyplot as plt
import h5py
import radmc



def make_config(model=1):
    if model == 1:
        config = radmc.StructuredJetModel.Config()
        config.table_resolution_radius = 1024
        config.table_resolution_theta = 256
        config.outermost_radius = 1e9
        config.jet_opening_angle = 0.3
        config.jet_polar_boundary = np.pi / 2
        config.jet_structure_exponent = 0.0
        config.specific_wind_power = 15
        config.luminosity_per_steradian = 1e49
        config.inner_radius_cm = 1e6
        config.leptons_per_baryon = 1.0
        config.photons_per_baryon = 1.0
        return config

    if model == 2:
        config = radmc.StructuredJetModel.Config()
        config.table_resolution_radius = 1024
        config.table_resolution_theta = 256
        config.outermost_radius = 1e9
        config.jet_opening_angle = 0.3
        config.jet_polar_boundary = np.pi / 2
        config.jet_structure_exponent = 0.0
        config.specific_wind_power = 24
        config.luminosity_per_steradian = 1e50
        config.inner_radius_cm = 1e6
        config.leptons_per_baryon = 1.0
        config.photons_per_baryon = 1.0
        return config

    if model == 3:
        config = radmc.StructuredJetModel.Config()
        config.table_resolution_radius = 1024
        config.table_resolution_theta = 256
        config.outermost_radius = 1e9
        config.jet_opening_angle = 0.3
        config.jet_polar_boundary = np.pi / 2
        config.jet_structure_exponent = 0.0
        config.specific_wind_power = 37
        config.luminosity_per_steradian = 1e51
        config.inner_radius_cm = 1e6
        config.leptons_per_baryon = 1.0
        config.photons_per_baryon = 1.0
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



def wien_photon(kT, nu):
    return 1. / (2 * kT**3) * nu**2 * np.exp(-nu / kT)



def wien_energy(kT, nu):
    return 1. / (2 * kT**3) * nu**3 * np.exp(-nu / kT)



class StructuredJetModel(object):
    def __init__(self, filename, restart=False, model=1):

        if restart:
            self.file = h5py.File(filename, 'a')
            self.photon_number = len(self.file['time'])
            config = self.read_config_from_file()

        else:
            config = make_config(model=model)
            self.file = h5py.File(filename, 'w')
            self.photon_number = 0
            self.write_config_to_file(config)

            for key in ['time', 'radius', 'theta', 'energy', 'lag']:
                if key not in self.file:
                    self.file.create_dataset(key, (0,), maxshape=(None,))

            self.file.create_group("tracks")

        self.config = config
        self.model = radmc.StructuredJetModel(config)
        self.photon_cache = []
        self.tracks_cache = []
        report_model(self.model, config)

    def run_photon(self, record_track=False):
        r0 = self.config.inner_radius_cm
        q = self.model.sample_theta(0.5)
        r = self.model.approximate_photosphere(q) / r0 * 0.01
        p = self.model.generate_photon(r, q)
        n = 0

        track = []

        while (p.position.radius < self.config.outermost_radius * 0.05
            and p.position.theta < self.config.jet_polar_boundary):
            p = self.model.step_photon(p)
            n += 1

            if record_track:
                track.append(p)

        self.photon_cache.append(p)
        self.photon_number += 1

        if record_track:
            self.tracks_cache.append(track)

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

        for track in self.tracks_cache:
            track_name = "{:06d}".format(len(self.file["tracks"]))
            track_group = self.file["tracks"].create_group(track_name)
            track_group['t'] = [p.position.t for p in track]
            track_group['x'] = [p.position.x for p in track]
            track_group['y'] = [p.position.y for p in track]
            track_group['z'] = [p.position.z for p in track]
            track_group['energy'] = [p.momentum.t for p in track]

        self.photon_cache = []
        self.tracks_cache = []
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



class WindProperties(object):

    def __init__(self):
        config = make_config()
        self.config = config
        self.model = radmc.StructuredJetModel(config)
        report_model(self.model, config)

    def plot_wind_photospheres(self):
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
            ax1.text(rphot * 0.75, uphot * 1.25, r"$\delta t = {:.1e}s$".format(lag))

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r"radius, cm")
        ax1.set_ylabel(r"four velocity (c)")
        plt.show()

    def plot_wind_profile(self):
        fig = plt.figure(figsize=[9, 6])
        ax1 = fig.add_subplot(1, 1, 1)

        r = np.logspace(0, 6, 128)
        cm = self.config.inner_radius_cm

        theta = 0.0
        states = [self.model.sample_wind_spherical(ri, theta) for ri in r]

        ax1.plot(r * cm, [s.u for s in states], label='Radial four velocity ($c$)')
        ax1.plot(r * cm, [s.temperature() for s in states], label=r'Comoving photon temperature ($m_e c^2$)')
        ax1.plot(r * cm, r**(-1.0) * 20)
        ax1.plot(r * cm, r**(-2./3))

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel(r"radius, cm")
        ax1.legend(loc='best')
        plt.show()

    def plot_wind_structure(self):
        num_axes = 3
        fig = plt.figure(figsize=[6, 8])
        axes = [fig.add_subplot(num_axes, 1, n + 1) for n in range(num_axes)]

        cm = self.config.inner_radius_cm
        physics = radmc.PhysicsConstants()
        thetas = np.linspace(0.0, 0.66, 24)
        structure = []

        for theta in thetas:
            rphot = self.model.approximate_photosphere(theta)
            state = self.model.sample_wind_spherical(rphot / cm, theta)
            power = self.model.angular_luminosity(theta)
            lag = self.model.approximate_lag_time(theta)
            props = dict(
                temp=state.temperature() * state.g * physics.gram_to_MeV(physics.me),
                power=power,
                lag=lag)

            structure.append(props)

        for i in range(num_axes - 1):
            axes[i].xaxis.set_major_formatter(plt.NullFormatter())

        axes[0].plot(thetas, [s['power'] for s in structure])
        axes[1].plot(thetas, [s['temp'] for s in structure])
        axes[2].plot(thetas, [s['lag'] for s in structure])

        axes[0].set_yscale('log')
        axes[1].set_yscale('log')
        axes[2].set_yscale('log')

        axes[0].set_ylabel(r"$L(\theta)$ (erg/s)")
        axes[1].set_ylabel(r"$h \nu$ (MeV)")
        axes[2].set_ylabel(r"$\delta t$ (s)")
        axes[2].set_xlabel(r"$\theta_{\rm obs}$")
        plt.show()



class PhotonProperties(object):

    def __init__(self, filename):
        self.model = StructuredJetModel(filename, restart=True)
        self.file = self.model.file
        self.T = self.file['time'][...]
        self.R = self.file['radius'][...]
        self.Q = self.file['theta'][...]
        self.E = self.file['energy'][...]
        self.L = self.file['lag'][...]

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

    def load_photons(self):
        return self.T, self.R, self.Q, self.E, self.L

    def load_photons_for_arrival_times(self, t0, t1):
        T, R, Q, E, L = self.load_photons()
        I = np.where((t0 < L) * (L <= t1))
        return T[I], R[I], Q[I], E[I], L[I]

    def load_photons_for_viewing_angles(self, theta0, theta1):
        T, R, Q, E, L = self.load_photons()
        I = np.where((theta0 < Q) * (Q <= theta1))
        return T[I], R[I], Q[I], E[I], L[I]

    def plot_light_curve(self):
        fig = plt.figure(figsize=[6, 8])
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)

        T, R, Q, E, L = self.load_photons()

        Lhist, Lbins = np.histogram(L, bins='auto')
        Fhist, Fbins = np.histogram(L, bins=Lbins, weights=E) # energy flux

        Lhist[Lhist == 0] = 1
        t = 0.5 * (Lbins[1:] + Lbins[:1])
        ax1.plot(t, Fhist / Lhist)

        ax2.hist(L[E < 1.25 * E.mean()], histtype='step', bins=Lbins, normed=False, label=r'$10 - 50 keV$')
        ax2.hist(L[E > 1.25 * E.mean()], histtype='step', bins=Lbins, normed=False, label=r'$50 - 300 keV$')
        ax1.set_ylabel(r"mean photon energy")
        ax2.set_ylabel(r"photon flux")
        ax2.set_xlabel(r"$t - t_c$")
        ax1.set_xlim(0.0, 6.0)
        ax2.set_xlim(0.0, 6.0)
        ax2.legend(loc='best')
        ax1.xaxis.set_major_formatter(plt.NullFormatter())
        plt.show()

    def plot_time_stats(self):
        fig = plt.figure(figsize=[10, 8])
        ax1 = fig.add_subplot(1, 1, 1)

        for t in [0.0, 1.0, 2.0]:
            T, R, Q, E, L = self.load_photons_for_arrival_times(t, t + 1.0)

            kT = 1. / 3 * E.mean()
            bins = np.logspace(np.log10(min(E)), np.log10(max(E)), 48)
            ax1.hist(E, bins=bins, histtype='step', normed=True, color='k', lw=t + 0.5)
            ax1.plot(bins, wien_photon(kT, bins), c='k', ls='--', lw=0.5)

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        plt.show()

    def plot_angle_stats(self):
        fig = plt.figure(figsize=[10, 8])
        ax1 = fig.add_subplot(2, 2, 1)
        ax2 = fig.add_subplot(2, 2, 2)
        ax3 = fig.add_subplot(2, 2, 3)
        ax4 = fig.add_subplot(2, 2, 4)

        angles = [0.0, 3.0]#, 0.2, 0.4, 0.6]

        for i in range(len(angles) - 1):
            T, R, Q, E, L = self.load_photons_for_viewing_angles(angles[i], angles[i + 1])

            kT = 1. / 3 * E.mean()
            bins = np.logspace(np.log10(min(E)), np.log10(max(E)), 100)
            ax1.hist(E, bins=bins, histtype='step', normed=True)
            ax1.plot(bins, wien_photon(kT, bins))

            Q0, NQ = self.make_log_histogram(Q)

            hist, bins = np.histogram(L, bins='auto')
            ax2.loglog(Q0, NQ / (2 * np.pi * np.sin(Q0)))
            ax3.hist(L, histtype='step', bins=bins, log=False)
            ax4.scatter(L, Q, s=.1)

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylabel(r"$dN/dE$")
        ax2.set_ylabel(r"$dN/d\Omega$")
        ax3.set_xlim(0.0, 10.0)
        ax3.set_ylabel(r"photon flux")
        ax3.set_xlabel(r"$t - t_c$")
        ax4.set_xlabel(r"lag time")
        ax4.set_ylabel(r"observer angle")
        plt.show()

    def plot_tracks(self):
        fig = plt.figure(figsize=[10, 8])
        ax1 = fig.add_subplot(1, 1, 1)
        r0 = self.model.config.inner_radius_cm
        q0 = self.model.config.jet_polar_boundary

        for track in self.file['tracks']:
            x = r0 * self.file['tracks'][track]['x'][:]
            z = r0 * self.file['tracks'][track]['z'][:]
            ax1.plot(x, z, c='k', marker='o', ms=0.5, lw='0.1')

        q = np.linspace(-q0, q0, 128)
        r = self.model.model.approximate_photosphere(0.0)
        x = r * np.sin(q)
        z = r * np.cos(q)
        ax1.plot(x, z)
        plt.show()


if __name__ == "__main__":
    choices = ['run',
    'light_curve',
    'angle_stats',
    'time_stats',
    'wind_profile',
    'wind_photospheres',
    'wind_structure',
    'tracks']
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=choices)
    parser.add_argument("--photons", type=int, default=1000)
    parser.add_argument("--model", type=int, default=1)
    parser.add_argument("--restart", action='store_true')
    parser.add_argument("--file", default="structured_jet.h5")
    args = parser.parse_args()

    if args.command == 'run':
        radmc.seed(str(datetime.datetime.now()))
        model = StructuredJetModel(args.file, restart=args.restart, model=args.model)

        while model.photon_number < args.photons:

            try:
                model.run_photon(record_track=False)
            except RuntimeError as e:
                print("Bad photon:", e)

            if len(model.photon_cache) >= 10:
                model.purge_photons()

    elif args.command == 'light_curve': PhotonProperties(args.file).plot_light_curve()
    elif args.command == 'angle_stats': PhotonProperties(args.file).plot_angle_stats()
    elif args.command == 'time_stats': PhotonProperties(args.file).plot_time_stats()
    elif args.command == 'tracks': PhotonProperties(args.file).plot_tracks()
    elif args.command == 'wind_profile': WindProperties().plot_wind_profile()
    elif args.command == 'wind_photospheres': WindProperties().plot_wind_photospheres()
    elif args.command == 'wind_structure': WindProperties().plot_wind_structure()


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

