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

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL
from chemreac import BANDED
from chemreac.integrate import run


"""
Demo of chemical reaction diffusion system.
"""

# <geometric_diffusion.png>

def main(tend=10.0, N=25, nt=30):
    x = np.linspace(0.1, 1.0, N+1)
    y0 = (x[0]/2+x[1:])**2

    geoms = (FLAT, SPHERICAL, CYLINDRICAL)
    geom_name = {FLAT: 'Flat', SPHERICAL: 'Spherical', CYLINDRICAL: 'Cylindrical'}

    fig = plt.figure()
    res = []

    for G in geoms:
        sys = ReactionDiffusion(1, [], [], [], N=N, D=[0.02], x=x, geom=G)
        res.append(run(sys, y0, t0=0, tend=tend, nt=nt))

    for i, G in enumerate(geoms):
        tout, yout, info = res[i]
        ax = fig.add_subplot(2,3,G+1, projection='3d')

        # create supporting points in polar coordinates
        T,X = np.meshgrid(x[0]/2+x[1:], tout)
        ax.plot_surface(T, X, yout, rstride=1, cstride=1, cmap=cm.YlGnBu_r)
        #ax.set_zlim3d(0, 1)
        if G == FLAT:
            ax.set_xlabel('x / m')
        else:
            ax.set_xlabel('r / m')
        ax.set_ylabel('time / s')
        ax.set_zlabel(r'C / mol*m**-3')
        ax.set_title(geom_name[G])


    for i, G in enumerate(geoms):
        tout, yout, info = res[i]
        if i == 0:
            ax = fig.add_subplot(2,3,3+i+1)
            for j in range(3):
                tout, yout, info =res[j]
                if j == 0:
                    yprim = yout
                elif j == 1:
                    yprim = yout*(x[1:]**3-x[:-1]**3)
                else:
                    yprim = yout*(x[1:]**2-x[:-1]**2)
                ybis = np.sum(yprim, axis=1)
                ax.plot(tout, ybis, label=str(j))
            ax.legend(loc='best')
            ax.set_title('Mass conservation')
        else:
            yout = yout - res[0][1] # difference
            ax = fig.add_subplot(2,3,3+G+1, projection='3d')

            # create supporting points in polar coordinates
            T,X = np.meshgrid(x[0]/2+x[1:], tout)
            ax.plot_surface(T, X, yout, rstride=1, cstride=1, cmap=cm.YlGnBu_r)
            #ax.set_zlim3d(0, 1)
            if G == FLAT:
                ax.set_xlabel('x / m')
            else:
                ax.set_xlabel('r / m')
            ax.set_ylabel('time / s')
            ax.set_zlabel(r'C / mol*m**-3')
            ax.set_title(geom_name[G] + ' minus ' + geom_name[0])

    plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)