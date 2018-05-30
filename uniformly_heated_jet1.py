

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

        self.axu = fig.add_subplot(4, 1, 1)
        self.axw = fig.add_subplot(4, 1, 2)
        self.axn = fig.add_subplot(4, 1, 3)
        self.axT = fig.add_subplot(4, 1, 4)

        for ax in [self.axu, self.axw, self.axn, self.axT]:
            ax.set_xscale('log')
            ax.set_yscale('log')

        self.axT.set_xlabel(r'Radius relative to inner boundary, $r / r_0$')
        self.axu.set_ylabel(r'Four-velocity, $u = \gamma \beta$')
        self.axw.set_ylabel(r'Thermal enthalpy, $w_{\rm th}$')
        self.axn.set_ylabel(r'Free enthalpy, $w_f$')
        self.axT.set_ylabel(r'$\theta = k T / m_e c^2$')

        plt.setp(self.axu.get_xticklabels(), visible=False)
        plt.setp(self.axw.get_xticklabels(), visible=False)
        plt.setp(self.axn.get_xticklabels(), visible=False)
        # fig.suptitle("$n_\gamma / n_p = 10^4$, $\mu_f / \mu_t = 1/10$")
        fig.suptitle("Typical solutions for turbulent GRB jets")

    def configure(self, state):
        return state \
        .set_luminosity_per_steradian(1e51) \
        .set_inner_radius_cm(1e6) \
        .set_leptons_per_baryon(1.0) \
        .set_photons_per_baryon(1e4) \

    def plot_wind(self, free_power, heating_rate, label, ls='-', lw=2, c=None):
        import radmc, numpy as np

        wind = radmc.RelativisticWind()
        wind.set_specific_wind_power(100, free_power)
        wind.set_initial_four_velocity(1.0)
        wind.set_heating_rate(heating_rate)
        solution = [self.configure(s) for s in wind.integrate_range(1e15)]

        r = [s.r for s in solution]
        u = [s.u for s in solution]
        n = [s.n for s in solution]
        w = [s.w for s in solution]
        T = [s.temperature() for s in solution]

        #self.axT.plot(r,  3 * np.array(r)**-0.300, ls='--', c='orange')
        #self.axT.plot(r, 10 * np.array(r)**-0.666, ls='--', c='k')

        self.axu.plot(r, u, label=label, ls=ls, lw=lw, c=c)
        self.axw.plot(r, w, label=label, ls=ls, lw=lw, c=c)
        self.axn.plot(r, n, label=label, ls=ls, lw=lw, c=c)
        self.axT.plot(r, T, label=label, ls=ls, lw=lw, c=c)

        self.axw.set_ylim(1e-3, 1e2)
        self.axn.set_ylim(1e-4, 1e2)
        #self.axT.set_ylim(1e-7, 1e2)



class TurbulencePropertiesOverRadiusPlot(object):

    """A class that runs wind dissipative hydrodynamic wind solutions.
    """

    def __init__(self, fig):
        import matplotlib.pyplot as plt

        self.axv = fig.add_subplot(4, 1, 1)
        self.axl = fig.add_subplot(4, 1, 2)
        self.axt = fig.add_subplot(4, 1, 3)
        self.axq = fig.add_subplot(4, 1, 4)

        for ax in [self.axv, self.axl, self.axt, self.axq]:
            ax.set_xscale('log')
            ax.set_yscale('log')

        self.axv.set_ylabel(r'Turbulent velocity, $v_0 / c$')
        self.axl.set_ylabel(r'Reynolds number, $\rm{Re}$')
        self.axt.set_ylabel(r'Optical depth, $\tau$')
        self.axq.set_ylabel(r'$\rm{Re} / \tau^{4/3}$')
        self.axq.set_xlabel(r'Radius relative to inner boundary, $r / r_0$')
        fig.suptitle("Turbulence properties for different heating rates")

    def configure(self, state):
        return state \
        .set_luminosity_per_steradian(1e51) \
        .set_inner_radius_cm(1e6) \
        .set_leptons_per_baryon(1.0) \
        .set_photons_per_baryon(1e4) \

    def plot_wind(self, free_power, heating_rate, label, ls='-', lw=2, c=None):
        import radmc, numpy as np

        wind = radmc.RelativisticWind()
        wind.set_specific_wind_power(100, free_power)
        wind.set_initial_four_velocity(1.0)
        wind.set_heating_rate(heating_rate)
        solution = [self.configure(s) for s in wind.integrate_range(1e15)]

        def eddy_optical_depth(s):
            return s.causally_connected_scale() / s.thomson_mean_free_path_comoving()

        def eddy_reynolds_number(s, v0):
            physics = radmc.PhysicsConstants()
            l0 = s.causally_connected_scale()
            nu = s.radiation_viscosity()
            return l0 * v0 * physics.c / nu

        chi = 1.0
        r = np.array([s.r for s in solution])
        n = np.array([s.n for s in solution])
        u = np.array([s.u for s in solution])
        v = (n * chi * heating_rate * u / (1 + u**2)**0.5)**(1./3)
        Re  = np.array([eddy_reynolds_number(s, v0) for s, v0 in zip(solution, v)])
        tau = np.array([eddy_optical_depth(s) for s in solution])
        self.axv.plot(r, v, label=label, ls=ls, lw=lw, c=c)
        self.axl.plot(r, Re , label=label, ls=ls, lw=lw, c=c)
        self.axt.plot(r, tau, label=label, ls=ls, lw=lw, c=c)
        self.axq.plot(r, Re / tau**(4./3), label=label, ls=ls, lw=lw, c=c)

        self.axl.axhline(1.0, linestyle='--', lw=1, alpha=0.1, c='k')
        self.axt.axhline(1.0, linestyle='--', lw=1, alpha=0.1, c='k')
        self.axq.axhline(1.0, linestyle='--', lw=1, alpha=0.1, c='k')



class DissipationRates(metaclass=Command):

    def __call__(self, **kwargs):
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=[6, 10])
        plot = DissipativeWindPlot(fig)
        plot.plot_wind(0., 0.0,  r'$\zeta=0$ (Ideal)', ls='--', lw=1)
        plot.plot_wind(50, 1e-8, r'$\zeta=10^{-6}$', ls='-.', lw=1)
        plot.plot_wind(50, 1e-3, r'$\zeta=10^{-3}$', lw=4./3)
        plot.plot_wind(50, 1e-2, r'$\zeta=10^{-2}$', lw=16./9)
        plot.plot_wind(50, 1e-1, r'$\zeta=10^{-1}$', lw=64./27)
        plot.plot_wind(50, 1e+0, r'$\zeta=1.0$', lw=256./81)
        plot.axu.legend(loc='best', ncol=2)
        fig.subplots_adjust(top=0.94, bottom=0.06)

        if kwargs['hardcopy']:
            plt.savefig("DissipativeJetSolutions.pdf")
        else:
            plt.show()



class TurbulenceProperties(metaclass=Command):

    def __call__(self, **kwargs):
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=[6, 10])
        plot = TurbulencePropertiesOverRadiusPlot(fig)
        plot.plot_wind(0., 0.0, None)
        plot.plot_wind(50, 1e-8, r'$\zeta=10^{-6}$', ls='-.', lw=1)
        plot.plot_wind(50, 1e-3, r'$\zeta=10^{-3}$', lw=4./3)
        plot.plot_wind(50, 1e-2, r'$\zeta=10^{-2}$', lw=16./9)
        plot.plot_wind(50, 1e-1, r'$\zeta=10^{-1}$', lw=64./27)
        plot.plot_wind(50, 1e+0, r'$\zeta=1.0$', lw=256./81)
        plot.axv.legend(loc='lower left', ncol=2)
        fig.subplots_adjust(top=0.94, bottom=0.06)

        if kwargs['hardcopy']:
            plt.savefig("TurbulenceProperties.pdf")            
        else:
            plt.show()



class HeatedJetModel(object):

    def __init__(self, filename, restart=False, heating_rate=0.0):
        import math, h5py, radmc

        if restart:
            self.file = h5py.File(filename, 'a')
            self.photon_number = len(self.file['time'])
            config = self.read_config_from_file()

        else:
            config = radmc.StructuredJetModel.Config()
            config.approximate_electron_energies_as_delta = True
            config.table_resolution_radius = 256
            config.table_resolution_theta = 1
            config.outermost_radius = 1e7
            config.jet_structure_exponent = 0
            config.specific_wind_power = 100
            config.luminosity_per_steradian = 1e52 / 4 / math.pi
            config.heating_rate = heating_rate * config.specific_wind_power
            config.inner_radius_cm = 1e7
            config.leptons_per_baryon = 1.0
            config.photons_per_baryon = 1e5

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
        self.report_model(self.model, config)

    def run_photon(self, record_track=False):
        rp = self.model.approximate_photosphere(0.0)
        r0 = self.config.inner_radius_cm
        q = self.model.sample_theta(1.0)
        r = rp / r0 * 0.00420
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
        import numpy as np

        fig = plt.figure(figsize=[7.35, 3.0])
        ax = fig.add_subplot(1, 1, 1)

        files = sorted(file.split(','), key=lambda f: HeatedJetModel(f, restart=True).file['config']['heating_rate'].value)
        colors = plt.get_cmap('plasma')(np.linspace(.2, .85, len(files) + 1)[:-1])

        for f, color in zip(files, colors):

            if f == 'xi1.h5':
                bins, kT = self.run(f, ax, color)
            else:
                self.run(f, ax, color)

        ax.plot(bins, [3e3 * self.wien_photon_spectrum(kT, nu) for nu in bins],
            ls=':',
            lw=2,
            c=color,
            label=r'Wien')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r"$h \nu_{\rm obs} \ ({\rm MeV})$")
        ax.set_ylabel(r"$dN / d\log\nu_{\rm obs}$")
        ax.set_xlim(1e-5, 3e3)
        ax.set_ylim(4, 2e3)
        ax.legend(loc='upper right', ncol=1)
        fig.subplots_adjust(left=0.08, bottom=0.16, right=0.95, top=0.95, wspace=0.4, hspace=0.2)

        if kwargs['hardcopy']:
            plt.savefig("PhotonSpectrum.pdf")
        else:
            plt.show()

    def run(self, file, ax, color):
        import matplotlib.pyplot as plt
        import numpy as np

        model = HeatedJetModel(file, restart=True)
        E = model.file['energy'][:] * 0.511
        w = np.ones_like(E) * 20000. / 300000.
        xi = model.file['config']['heating_rate'].value / 100
        label = r"$\tilde \xi={}$".format(self.xi_to_string(xi))

        bins = np.logspace(np.log10(min(E)), np.log10(max(E)), 100)
        ax.hist(E, bins=bins, weights=w, histtype='step', label=label, color=color, lw=0.5)
        ax.hist(E, bins=bins, weights=w, histtype='stepfilled', color=color, alpha=0.1)
        return bins, E.mean() / 3

    def wien_photon_spectrum(self, kT, nu):
        import math
        return 1. / (2 * kT**3) * nu**2 * math.exp(-nu / kT)

    def xi_to_string(self, xi):
        import numpy as np
        if xi == 0.0:
            return r"0"
        if xi == 1.0:
            return r"1"
        return r"10^{{{:.0f}}}".format(np.log10(xi))



if __name__ == "__main__":
    from matplotlib import rc
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('command', choices=Command.get_commands())
    parser.add_argument("--photons", type=int, default=1000)
    parser.add_argument("--restart", action='store_true')
    parser.add_argument("--heating_rate", type=float, default=0.0)
    parser.add_argument("--hardcopy", action='store_true')
    parser.add_argument("--file", type=str, default='heated_jet.h5')
    args = parser.parse_args()

    rc('text', usetex=args.hardcopy)
    rc('font', **{'size':9})
    rc('xtick', labelsize=9)
    rc('ytick', labelsize=9)
    rc('axes', labelsize=9)

    cmd = Command.get_command_class(args.command)()
    cmd(**vars(args))
