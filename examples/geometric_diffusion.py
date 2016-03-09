#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import argh
import numpy as np

from chemreac import ReactionDiffusion, Geom_names
from chemreac.integrate import run
from chemreac.util.plotting import plot_C_vs_t_and_x, save_and_or_show_plot

"""
Demo of diffusion.
"""

# <geometric_diffusion.png>


def main(tend=10.0, N=25, nt=30, nstencil=3, linterpol=False,
         rinterpol=False, num_jacobian=False, plot=False,
         savefig='None', verbose=False):
    x = np.linspace(0.1, 1.0, N+1)

    def f(x):
        return 2*x**2/(x**4+1)  # f(0)=0, f(1)=1, f'(0)=0, f'(1)=0
    y0 = f(x[1:])+x[0]  # (x[0]/2+x[1:])**2

    t0 = 1e-10
    tout = np.linspace(t0, tend, nt)

    res, systems = [], []
    for g in 'fcs':
        rd = ReactionDiffusion(1, [], [], [], N=N, D=[0.02], x=x,
                               geom=g, nstencil=nstencil, lrefl=not linterpol,
                               rrefl=not rinterpol)
        integr = run(rd, y0, tout, with_jacobian=(not num_jacobian))
        res.append(integr.yout)
        systems.append(rd)

    if plot:
        # matplotlib
        from mpl_toolkits.mplot3d import Axes3D
        assert Axes3D  # silence pyflakes
        from matplotlib import cm
        from matplotlib import pyplot as plt

        fig = plt.figure()

        # Plot spatio-temporal conc. evolution
        for i, g in enumerate('fcs'):
            yout = res[i]
            ax = fig.add_subplot(2, 3, 'fcs'.index(g)+1, projection='3d')

            plot_C_vs_t_and_x(rd, tout, yout, 0, ax,
                              rstride=1, cstride=1, cmap=cm.gist_earth)
            ax.set_title(Geom_names[g])

        # Plot mass conservation
        ax = fig.add_subplot(2, 3, 4)
        for j in range(3):
            ax.plot(tout, np.apply_along_axis(
                systems[j].integrated_conc, 1, res[j][:, :, 0]), label=str(j))
        ax.legend(loc='best')
        ax.set_title('Mass conservation')

        # Plot difference from flat evolution (not too informative..)
        for i, g in enumerate('fcs'):
            yout = res[i]  # only one specie
            if i != 0:
                yout = yout - res[0]  # difference (1 specie)
                ax = fig.add_subplot(2, 3, 3+'fcs'.index(g)+1, projection='3d')

                plot_C_vs_t_and_x(rd, tout, yout, 0, ax,
                                  rstride=1, cstride=1, cmap=cm.gist_earth)
                ax.set_title(Geom_names[g] + ' minus ' + Geom_names[0])

        save_and_or_show_plot(savefig)


if __name__ == '__main__':
    argh.dispatch_command(main)
