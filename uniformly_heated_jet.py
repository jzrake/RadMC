


class Command(type):
    commands = dict()

    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)
        meta.commands[meta.convert(name)] = cls
        return cls

    def get_commands():
        return Command.commands.keys()

    def get_command_class(name):
        return Command.commands[name]

    def convert(name):
        import re
        s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()



class DissipativeWindPlot(object):

    """A class that runs wind dissipative hydrodynamic wind solutions.
    """

    def __init__(self, fig):
        import matplotlib.pyplot as plt

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
        fig.suptitle("$n_\gamma / n_p = 10^4$, $\mu_f / \mu_t = 1/10$")

    def configure(self, state):
        return state.set_photons_per_baryon(1e4)

    def plot_wind(self, free_enthalpy, heating_rate, label):
        import radmc

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
        if True or heating_rate > 0:
            self.axn.plot(r, n, label=label)



class DissipationRates(metaclass=Command):

    def __call__(self, **kwargs):
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=[6, 10])
        plot = DissipativeWindPlot(fig)
        plot.plot_wind(0.00, 0.0, 'ideal')
        plot.plot_wind(1, 0.01, r'$\zeta=0.01$')
        plot.plot_wind(1, 0.1, r'$\zeta=0.1$')
        plot.plot_wind(1, 0.2, r'$\zeta=0.2$')
        plot.plot_wind(1, 0.4, r'$\zeta=0.4$')
        plot.plot_wind(1, 0.8, r'$\zeta=0.8$')
        plot.axu.legend(loc='best')

        plt.show()



class HeatedJetModel(object):

    def __init__(self, filename, restart=False, heating_rate=0.0):
        import h5py, radmc

        if restart:
            self.file = h5py.File(filename, 'a')
            self.photon_number = len(self.file['time'])
            config = self.read_config_from_file()

        else:
            config = radmc.StructuredJetModel.Config()
            config.approximate_electron_energies_as_delta = True
            config.table_resolution_radius = 256
            config.table_resolution_theta = 1
            config.outermost_radius = 1e9
            config.jet_structure_exponent = 0
            config.specific_wind_power = 100
            config.luminosity_per_steradian = 1e49
            config.heating_rate = heating_rate
            config.initial_free_enthalpy = 10.0 if heating_rate > 0.0 else 0.0
            config.inner_radius_cm = 1e6
            config.leptons_per_baryon = 1.0
            config.photons_per_baryon = 1e5

            self.file = h5py.File(filename, 'w')
            self.photon_number = 0
            self.write_config_to_file(config)

            for attr in ['initial_free_enthalpy', 'heating_rate']:
                print("{} ... {}".format(attr, getattr(config, attr)))

            for key in ['time', 'radius', 'theta', 'energy', 'lag']:
                if key not in self.file:
                    self.file.create_dataset(key, (0,), maxshape=(None,))

            self.file.create_group("tracks")

        self.config = config
        self.model = radmc.StructuredJetModel(config)
        self.photon_cache = []
        self.tracks_cache = []
        self.report_model(self.model, config)

    def run_photon(self, record_track=False):
        rp = self.model.approximate_photosphere(0.0)
        r0 = self.config.inner_radius_cm
        q = self.model.sample_theta(1.0)
        r = rp / r0 * 0.01
        p = self.model.generate_photon(r, q)
        n = 0

        track = []

        while p.position.radius < rp / r0 * 2:
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
        import radmc

        config_group = self.file['config']
        config = radmc.StructuredJetModel.Config()

        for key in dir(config):
            if key.startswith('__'): continue
            try:
                setattr(config, key, config_group[key].value)
            except KeyError as e:
                print(e)
        return config

    def report_model(self, model, config):
        import radmc

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



class PropagatePhotons(metaclass=Command):

    def __call__(self, restart='', heating_rate=0.0, **kwargs):
        import matplotlib.pyplot as plt
        import numpy as np

        model = HeatedJetModel(kwargs['file'], restart=restart, heating_rate=heating_rate)

        for i in range(kwargs['photons']):
            model.run_photon(record_track=False)

            if len(model.photon_cache) >= 100:
                model.purge_photons()



class HistogramPhotons(object, metaclass=Command):

    def __call__(self, file='', **kwargs):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        for f in file.split(','):
            bins, T = self.run(f, ax)

        ax.plot(bins, [1e2 * e * self.wien_photon_spectrum(T, e) for e in bins], label='Wien')
        ax.plot(bins, 5e4 * bins**1.4, label=r'$E^{1.4}$', lw=0.5, ls='--', c='k')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r"$h \nu$")
        ax.set_ylabel(r"$dN / d\log\nu$")
        ax.set_ylim(1e-3, 1e2)
        ax.legend(loc='best')
        plt.show()

    def run(self, file, ax):
        import matplotlib.pyplot as plt
        import numpy as np

        model = HeatedJetModel(file, restart=True)
        E = model.file['energy'][:]

        T = 1. / 3 * E.mean()
        bins = np.logspace(np.log10(min(E)), np.log10(max(E)), 100)
        ax.hist(E, bins=bins, density=True, histtype='step', weights=E, label=file)
        # ax.plot(bins, [1e2 * e * self.wien_photon_spectrum(T, e) for e in bins], label='Wien')

        return bins, T

    def wien_photon_spectrum(self, kT, nu):
        import math
        return 1. / (2 * kT**3) * nu**2 * math.exp(-nu / kT)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('command', choices=Command.get_commands())
    parser.add_argument("--photons", type=int, default=1000)
    parser.add_argument("--restart", action='store_true')
    parser.add_argument("--heating_rate", type=float, default=0.0)
    parser.add_argument("--file", type=str, default='heated_jet.h5')
    args = parser.parse_args()

    cmd = Command.get_command_class(args.command)()
    cmd(**vars(args))
