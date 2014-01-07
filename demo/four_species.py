#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argh
import numpy as np

from chemreac.serialization import load
from chemreac.cpp_chem_wrapper import \
    PyReactionDiffusion as ReactionDiffusion

from chemreac import DENSE, BANDED, SPARSE
from chemreac.integrate import run

import matplotlib.cm
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

"""
Demo of chemical reaction diffusion system.
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0

def coloured_spy(A, cmap_name='gray'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.imshow(A, cmap=matplotlib.cm.get_cmap(cmap_name),
               interpolation='none')
    ax = plt.gca()
    ya = ax.get_yaxis()
    ya.set_major_locator(MaxNLocator(integer=True))
    xa = ax.get_xaxis()
    xa.set_major_locator(MaxNLocator(integer=True))
    plt.colorbar()


def main(tend=10.0, N=1, nt=50, plot=True, spy=False, mode=None):
    sys = load('four_species.json', N=N)

    y0 = [1.3, 1e-4, 0.7, 1e-4]
    y0 = np.array(y0*sys.N)
    #y0 = np.array([[float(x)]*4 for x in range(1,N+1)]).flatten()

    t0 = 0.0

    if mode == None:
        if sys.N == 1:
            mode = DENSE
        elif sys.N > 1:
            mode = BANDED
    else:
        mode = int(mode)

    if spy:
        fout = np.empty(sys.n*sys.N)
        sys.f(t0, y0, fout)
        print(fout)
        if mode == DENSE:
            jout = np.zeros((sys.n*sys.N, sys.n*sys.N), order='F')
            sys.dense_jac_cmaj(t0, y0, jout)
            coloured_spy(np.log(jout))
        elif mode == BANDED:
            jout = np.zeros((sys.n*2+1, sys.n*sys.N), order='F')
            sys.banded_packed_jac_cmaj(t0, y0, jout)
            coloured_spy(np.log(jout))
        plt.show()

    else:

        tout, yout, info = run(sys, y0, t0, tend, nt)
        if plot:
            for i,l in enumerate('ABCD'):
                plt.plot(tout, yout[:,i], label=l)
            plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
