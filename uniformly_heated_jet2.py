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



def xi_to_string(xi):
    import numpy as np
    if xi == 0.0:
        return r"0"
    if xi == 1.0:
        return r"1"
    return r"10^{{{:.0f}}}".format(np.log10(xi))



def wave_temperature(s, beta0):
    """
    Compute the turbulent wave temperature.
    """
    import numpy as np

    if beta0 == 0.0:
        return 0.0

    tau = s.jet_optical_depth() 
    Re = beta0 * tau / s.w
    ell0 = s.causally_connected_scale()
    ellst = s.thomson_mean_free_path_comoving()
    ellnu = ell0 * Re**(-3. / 4)
    velnu = beta0 * (ellnu / ell0)**(1. / 3)

    theta_w_visc = (1. / 3) * velnu**2 * (ellst / ellnu)**2
    theta_w_turb = (1. / 3) * beta0**2 * (ellst / ell0)**(2. / 3)

    x = np.log(ellst / ellnu)
    y = (1 + np.tanh(x)) / 2
    # y = 1 if x > 0 else 0

    return theta_w_turb * y + theta_w_visc * (1 - y)



def uniformly_heated_jet_solution(
    eta0=100.0,
    heating_rate=0.0,
    inner_radius=1e7,
    outer_radius=1e14,
    photons_per_baryon=1e5,
    just_return_labels=False):
    """
    Evolve a uniformly heated jet model for the given parameters. Return a
    dictionary of 1D arrays containing thermodynamic and kinematic jet
    properties.
    """
    import math, radmc

    if just_return_labels:
        return dict(
            r=r'$r \, \, (\rm{cm})$',
            u=r'$u$',
            w=r'$w$',
            y=r'$y$',
            tau=r'$\tau$',
            theta_e=r'$\theta_e$',
            theta_g=r'$\theta_\gamma$',
            theta_w=r'$\theta_w$',
            d_theta=r'$\Delta \theta$',
            ell0=r"$\ell_0$",
            nu=r"$\nu$",
            Re=r"${\rm Re}$",
            Retau43=r"${\rm Re} / \tau^{4/3}$",
            TwTe=r'$\theta_w / \theta_{\rm eff}$',
            TwTt=r'$\theta_w / \theta_e$')

    def configure_state(state, wind):
        return state \
        .set_luminosity_per_steradian(1e52 / 4 / math.pi) \
        .set_inner_radius_cm(inner_radius) \
        .set_leptons_per_baryon(1.0) \
        .set_photons_per_baryon(photons_per_baryon)

    wind = radmc.RelativisticWind()
    wind.set_specific_wind_power(eta0)
    wind.set_initial_four_velocity(10.0)
    wind.set_heating_rate(heating_rate * eta0) # Note: this is appropriate for the xi-tilde prescription
    solution = [configure_state(s, wind) for s in wind.integrate_range(outer_radius / inner_radius)]
    beta0 = heating_rate**(1. / 3)

    r = [s.radius() for s in solution]
    u = [s.u for s in solution]
    w = [s.w for s in solution]
    y = [s.compton_parameter() for s in solution]
    tau = [s.jet_optical_depth() for s in solution]
    theta_e = [s.electron_temperature() for s in solution] # effective electron temp
    theta_g = [s.photon_temperature() for s in solution]
    theta_w = [wave_temperature(s, beta0) for s in solution]
    theta_t = [s.electron_temperature() - wave_temperature(s, beta0) for s in solution] # true electron temp
    d_theta = [s.delta_theta() for s in solution]
    ell0    = [s.causally_connected_scale() for s in solution]
    nu      = [s.radiation_viscosity() for s in solution]
    Re      = [beta0 * s.jet_optical_depth() / s.w for s in solution]
    Retau43 = [beta0 * s.jet_optical_depth()**(-1. / 3) / s.w for s in solution]
    TwTe    = [tw / te for te, tw in zip(theta_e, theta_w)]
    TwTt    = [tw / tt for tt, tw in zip(theta_t, theta_w)]

    return dict(r=r, u=u, w=w, y=y,
        tau=tau,
        theta_e=theta_e, # effective electron temperature (wave plus thermal)
        theta_g=theta_g, # Compton temperature
        theta_w=theta_w, # wave temperature
        theta_t=theta_t, # thermal electron temperature
        d_theta=d_theta,
        ell0=ell0,
        nu=nu,
        Re=Re,
        Retau43=Retau43,
        TwTe=TwTe,
        TwTt=TwTt)



class ThermodynamicEvolution(object, metaclass=Command):

    def __call__(self, file='', **kwargs):
        import matplotlib.pyplot as plt
        import numpy as np

        L = uniformly_heated_jet_solution(just_return_labels=True)

        xis = [0.0, 1e-3, 1e-2, 1e-1, 1.0]
        results = [uniformly_heated_jet_solution(1e2, heating_rate=xi) for xi in xis]

        series_names = ['u', 'w', 'y', 'theta_g', 'theta_e', 'tau']
        series_datas = [[s[n] for n in series_names] for s in results]
        theta_es = [series_datas[n][4] for n in range(len(xis))]
        taus     = [series_datas[n][5] for n in range(len(xis))]
        xdatas   = [s['r'] for s in results]
        colors   = plt.get_cmap('plasma')(np.linspace(.4, .85, len(xis)))

        fig = plt.figure(figsize=[7.35, 5.50])
        axes = [fig.add_subplot(2, 2, n) for n in [1, 2, 3, 4]]

        for xdata, series_data, xi, color in zip(xdatas, series_datas, xis, colors):
            for ax, name, ydata in zip(axes, series_names, series_data):

                if max(ydata) < 1e-10: continue # to skip zero y data

                label = r"$\tilde \xi={}$".format(xi_to_string(xi))
                ax.plot(xdata, ydata, label=label, color=color)

        for n in range(len(xis)):
            axes[3].plot(xdatas[n], theta_es[n], ls='-', lw=3, alpha=0.3, color=colors[n])

        photosphere_indexes = [np.where((np.array(taus[n]) < 1.0))[0][0] for n in range(len(xis))]

        for m, ax in enumerate(axes):
            for n, i in enumerate(photosphere_indexes):

                if ax is axes[2] and n == 0: continue

                ax.scatter(xdatas[n][i], series_datas[n][m][i], s=15, color=colors[n])

        for ax, name in zip(axes, series_names):
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylabel(L[name])
            if ax is axes[2] or ax is axes[3]:
                ax.set_xlabel(L['r'])

        axes[0].legend(loc='best', ncol=2)
        axes[2].set_ylim(1e-3, 1e1)
        axes[2].axhline(1.0, ls=':', lw=1, c='grey')
        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.3, hspace=0.2)

        if kwargs['hardcopy']:
            plt.savefig('ThermodynamicEvolution.pdf')
        else:
            plt.show()



class ThermodynamicEvolution2(object, metaclass=Command):

    def __call__(self, file='', **kwargs):
        import matplotlib.pyplot as plt
        import numpy as np

        L = uniformly_heated_jet_solution(just_return_labels=True)
        L["u"] = r"$\Gamma$"

        xis = [0.0, 1e-3, 1e-2, 1e-1, 1.0]
        results = [uniformly_heated_jet_solution(1e2, heating_rate=xi) for xi in xis]

        series_names = ['u', 'w', 'y', 'theta_g', 'theta_e', 'tau']
        colors   = plt.get_cmap('plasma')(np.linspace(.4, .85, len(xis)))

        fig = plt.figure(figsize=[7.35, 5.50])
        axes = [fig.add_subplot(2, 2, n) for n in [1, 2, 3, 4]]

        for n in range(len(xis)):
            label = r"$\tilde \xi={}$".format(xi_to_string(xis[n]))

            r  = np.array(results[n]['r'])
            u  = np.array(results[n]['u'])
            w  = np.array(results[n]['w'])
            y  = np.array(results[n]['y'])
            Tg = np.array(results[n]['theta_g'])
            Te = np.array(results[n]['theta_e'])

            if n == 1:
                ax3labelTg = r"Compton temp. $\theta_c$"
                ax3labelTe = r"Effective temp. $\theta_{\rm eff}$"
            else:
                ax3labelTg = None
                ax3labelTe = None

            # indexes below photosphere
            ithick = np.where(np.array(results[n]['tau']) > 1.0)

            axes[0].plot(r[ithick], u[ithick], c=colors[n], lw=2, label=label)
            axes[1].plot(r[ithick], w[ithick], c=colors[n], lw=2)
            axes[2].plot(r[ithick], y[ithick], c=colors[n], lw=2)
            axes[3].plot(r[ithick], Tg[ithick], c=colors[n], label=ax3labelTg, lw=2)
            axes[3].plot(r[ithick], Te[ithick], c=colors[n], label=ax3labelTe, lw=1, ls='--')

            axes[0].scatter(r[ithick][-1:], u[ithick][-1:], color=colors[n])
            axes[1].scatter(r[ithick][-1:], w[ithick][-1:], color=colors[n])
            axes[2].scatter(r[ithick][-1:], y[ithick][-1:], color=colors[n])
            axes[3].scatter(r[ithick][-1:], Tg[ithick][-1:], color=colors[n])

            if n == 1:
                m = -30
                axes[3].vlines(r[ithick][m:m+1],
                    ymin=Tg[ithick][m],
                    ymax=Te[ithick][m],
                    lw=1,
                    color='k',
                    linestyle='-',
                    zorder=100)

                axes[3].hlines(Te[ithick][m],
                    xmin=r[ithick][m] / 1.3,
                    xmax=r[ithick][m] * 1.3,
                    lw=1,
                    color='k',
                    linestyle='-',
                    zorder=100)

                axes[3].hlines(Tg[ithick][m],
                    xmin=r[ithick][m] / 1.3,
                    xmax=r[ithick][m] * 1.3,
                    lw=1,
                    color='k',
                    linestyle='-',
                    zorder=100)

                axes[3].text(r[ithick][m:m+1] * 1.5, 1e-3, r"$\Delta \theta$")

        for ax, name in zip(axes, series_names):
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylabel(L[name])
            if ax is axes[2] or ax is axes[3]:
                ax.set_xlabel(L['r'])


        axes[0].legend(loc='best', ncol=2)
        axes[2].set_ylim(5e-4, 1e1)
        axes[2].axhline(1.0, ls=':', lw=1, c='grey')

        axes[3].set_ylabel(r'$\theta_c$ (and $\theta_{\rm eff}$)')
        axes[3].legend(loc='best')

        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95, wspace=0.3, hspace=0.2)

        if kwargs['hardcopy']:
            plt.savefig('ThermodynamicEvolution.pdf')
        else:
            plt.show()



class TurbulenceEvolution2(object, metaclass=Command):

    def __call__(self, file='', **kwargs):
        import matplotlib.pyplot as plt
        import numpy as np

        L = uniformly_heated_jet_solution(just_return_labels=True)
        series_names = ['TwTt']
        # xis = [1e-3, 1e-2, 1e-1, 1.0]
        xis = [1e-3]
        results = [uniformly_heated_jet_solution(1e2, heating_rate=xi) for xi in xis]

        fig = plt.figure(figsize=[3.55, 3.55])
        ax = fig.add_subplot(1, 1, 1)
        colors = plt.get_cmap('plasma')(np.linspace(.2, .85, len(xis) + 1)[1:])

        for n in range(len(xis)):
            label = r"$\tilde \xi={}$".format(xi_to_string(xis[n]))

            r       = np.array(results[n]['r'])
            Retau43 = np.array(results[n]['Retau43'])
            tau     = np.array(results[n]['tau'])

            Tg      = np.array(results[n]['theta_g'])
            Te      = np.array(results[n]['theta_e'])
            Tw      = np.array(results[n]['theta_w']) / 1 # <-- NOTE: fudged wave temp
            Tt      = Te - Tw

            # ax.plot(r, Tw / (Tt - Tg), color=colors[n], label=label)
            ax.plot(r, Tw, lw=3, label="wave")
            ax.plot(r, Tt - Tg, lw=1, label="thermal minus Comp")
            ax.plot(r, Tg, lw=1, label="Compton")
            ax.plot(r, Te, lw=1, label="effective")

        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_ylabel(r"$\theta_w / \theta_e$")
        ax.set_xlabel(L['r'])
        ax.set_xlim(2e+7, 5e13)
        # ax.set_ylim(2e-3, 6e2)
        ax.legend(loc='best')
        fig.subplots_adjust(left=0.16, bottom=0.12, right=0.95, top=0.95, wspace=0.4, hspace=0.2)

        if kwargs['hardcopy']:
            plt.savefig('TurbulenceEvolution.pdf')
        else:
            plt.show()



class TurbulenceEvolution(object, metaclass=Command):

    def __call__(self, file='', **kwargs):
        import matplotlib.pyplot as plt
        import numpy as np

        L = uniformly_heated_jet_solution(just_return_labels=True)
        series_names = ['TwTe']
        xis = [1e-3, 1e-2, 1e-1, 1.0]
        results = [uniformly_heated_jet_solution(1e2, heating_rate=xi) for xi in xis]

        fig = plt.figure(figsize=[3.55, 3.55])
        ax = fig.add_subplot(1, 1, 1)
        colors = plt.get_cmap('plasma')(np.linspace(.2, .85, len(xis) + 1)[1:])

        for n in range(len(xis)):
            label = r"$\tilde \xi={}$".format(xi_to_string(xis[n]))

            r       = np.array(results[n]['r'])
            Retau43 = np.array(results[n]['Retau43'])
            tau     = np.array(results[n]['tau'])
            TwTe    = np.array(results[n]['TwTe'])

            # indexes in viscous limit
            ivisc = np.where(Retau43 < 1.0)

            # indexes where waves are hotter
            iwave = np.where(TwTe > 1.0)

            # indexes in turbulent limit below photosphere
            iturb = np.where((Retau43 > 1.0) * (tau > 1.0))

            # indexes above photosphere
            ithin = np.where(tau < 1.0)

            ax.plot(r[ivisc], TwTe[ivisc], c=colors[n], lw=3, label=label)
            ax.plot(r[iturb][4:], TwTe[iturb][4:], c=colors[n], lw=1)
            ax.plot(r[ithin][6:], TwTe[ithin][6:], c=colors[n], lw=1, alpha=1.0, ls='--')
            ax.fill_between(r[iwave], 1.0, TwTe[iwave], color=colors[n], alpha=0.25)
            ax.scatter(r[ithin][0], TwTe[ithin][0], color=colors[n])

        turb_label  = r"$\rightarrow$ $\ell_\star > \ell_\nu$"
        shear_label = r"$\ell_\star < \ell_\nu$ $\leftarrow$"

        ax.vlines([r[iturb][0]], 1, 25, zorder=20, lw=1)
        ax.text(r[iturb][0] * 1.2, 20, turb_label, horizontalalignment='left')
        ax.text(r[iturb][0] / 1.3, 20, shear_label, horizontalalignment='right')

        ax.annotate(
            'Turbulent Comp.',
            xy=(2e10, 1.5),
            xytext=(2e11, 5),
            horizontalalignment='left',
            arrowprops=dict(arrowstyle='->'))

        ax.annotate(
            'Shear Comp.',
            xy=(2e9, 1.5),
            xytext=(1e9, 5),
            horizontalalignment='right',
            arrowprops=dict(arrowstyle='->'))

        ax.text(2e11, 1e-2,
            'Thermal Comp.',
            horizontalalignment='left')

        ax.axhline(1.0, lw=1.0, ls=':', c='k')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylabel(L['TwTe'])
        ax.set_xlabel(L['r'])
        ax.set_xlim(2e+7, 5e13)
        ax.set_ylim(2e-3, 6e2)
        ax.legend(loc='upper right')
        fig.subplots_adjust(left=0.16, bottom=0.12, right=0.95, top=0.95, wspace=0.4, hspace=0.2)

        if kwargs['hardcopy']:
            plt.savefig('TurbulenceEvolution.pdf')
        else:
            plt.show()



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
