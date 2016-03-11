#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import os
from pprint import pprint

import argh
import numpy as np

from chemreac.serialization import load
from chemreac.integrate import run


"""
Demo of chemical reaction diffusion system.
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0


def main(tend=3.0, N=30, nt=30, plot=False, logy=False, logt=False,
         savefig='None', verbose=False):
    def mod1(x):
        return x/(x**2+1)

    # decay A->B is modulated with x
    rd = load(
        os.path.join(os.path.dirname(
            __file__), 'four_species.json'),
        N=N, modulation=[[mod1(x/3+0.1) for x in range(N)]],
        modulated_rxns=[1], logt=logt, logy=logy)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    # y0 = np.concatenate([y0/(i+1)*(0.25*i**2+1) for i in range(N)])
    y0 = np.concatenate([y0 for i in range(N)])

    t0 = 1e-10
    tout = np.linspace(t0, tend, nt)
    integr = run(rd, y0, tout)
    if verbose:
        pprint(integr.info)

    cx = rd.x[:-1]+np.diff(rd.x)/2  # center x
    if plot:
        # Using matplotlib for 3D plots:
        from mpl_toolkits.mplot3d import Axes3D
        assert Axes3D  # silence pyflakes
        from matplotlib import cm
        from matplotlib import pyplot as plt

        from chemreac.util.plotting import save_and_or_show_plot

        fig = plt.figure()

        for i, l in enumerate('ABCD'):
            ax = fig.add_subplot(3, 2, i+1, projection='3d')
            T, X = np.meshgrid(cx, tout)
            ax.plot_surface(T, X, integr.Cout[:, :, i], rstride=1, cstride=1,
                            cmap=cm.YlGnBu_r)
            ax.set_xlabel('x / m')
            ax.set_ylabel('time / s')
            ax.set_zlabel(r'C / mol*m**-3')
            ax.set_title(l)
            ax.legend(loc='best')

        ax = fig.add_subplot(3, 2, 5)
        print(np.array(rd.modulation).shape)
        ax.plot(cx, np.array(rd.modulation)[0, :])

        save_and_or_show_plot(savefig=savefig)

if __name__ == '__main__':
    argh.dispatch_command(main)
