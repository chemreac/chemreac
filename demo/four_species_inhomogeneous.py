#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import argh
import numpy as np
# matplotlib
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import cm
from matplotlib import pyplot as plt

from chemreac.serialization import load
from chemreac import DENSE, BANDED, SPARSE
from chemreac.integrate import run
from chemreac.util import coloured_spy


"""
Demo of chemical reaction diffusion system.
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0


def main(tend=3.0, N=30, nt=30, plot=False, mode=None, 
         logy=False, logt=False, show=False):
    mod1 = lambda x: x/(x**2+1)
    sys = load('four_species.json', N=N,
               bin_k_factor=[[mod1(x/3+0.1)] for x in range(N)], # decay A->B is modulated with x
               bin_k_factor_span=[1]
    )

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    y0 = np.concatenate([y0/(i+1)*(0.25*i**2+1) for i in range(N)])

    t0 = 1e-10
    tout = np.linspace(t0, tend, nt)
    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(sys, y, t)
    if logy: yout = np.exp(yout)

    cx = sys.x[:-1]+np.diff(sys.x)/2 # center x
    if plot:
        fig = plt.figure()

        for i,l in enumerate('ABCD'):
            ax = fig.add_subplot(3, 2, i+1, projection='3d')
            T, X = np.meshgrid(cx, tout)
            ax.plot_surface(T, X, yout[:,i::4], rstride=1, cstride=1, cmap=cm.YlGnBu_r)
            ax.set_xlabel('x / m')
            ax.set_ylabel('time / s')
            ax.set_zlabel(r'C / mol*m**-3')
            ax.set_title(l)
            ax.legend(loc='best')

        ax = fig.add_subplot(3, 2, 5)
        print(sys.bin_k_factor[:,0].shape)
        ax.plot(cx, sys.bin_k_factor[:,0])

        if show:
            plt.show()
        else:
            plt.savefig(__file__[:-2]+'png')

if __name__ == '__main__':
    argh.dispatch_command(main)
