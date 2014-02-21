#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac.serialization import load
from chemreac import DENSE, BANDED, SPARSE
from chemreac.integrate import run
from chemreac.util import coloured_spy


"""
Demo of chemical reaction diffusion system.
"""

# A      -> B          k1=0.05
# 2C + B -> D + B      k2=3.0


def main(tend=10.0, N=1, nt=500, plot=False, jac_spy=False, mode=None,
         logy=False, logt=False, show=False):

    sys = load('four_species.json', N=N, x=N, logy=logy, logt=logt)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    y0 = np.concatenate([y0/(i+1)*(0.25*i**2+1) for i in range(N)])
    t0 = 1e-10

    if mode == None:
        if sys.N == 1:
            mode = DENSE
        elif sys.N > 1:
            mode = BANDED
    else:
        mode = int(mode)

    if jac_spy:
        fout = np.empty(sys.n*sys.N)
        sys.f(t0, y0, fout)
        print(fout)
        if mode == DENSE:
            jout = np.zeros((sys.n*sys.N, sys.n*sys.N), order='F')
            sys.dense_jac_cmaj(t0, y0, jout)
            coloured_spy(np.log(jout))
        elif mode == BANDED:
            jout = np.zeros((sys.n*2+1, sys.n*sys.N), order='F') # note sys.n*3 needed in actual call
            sys.banded_packed_jac_cmaj(t0, y0, jout)
            coloured_spy(np.log(jout))
        plt.show()

    else:

        tout = np.linspace(t0, tend, nt)
        y = np.log(y0) if logy else y0
        t = np.log(tout) if logt else tout
        yout, info = run(sys, y, t)
        if logy: yout = np.exp(yout)
        if plot:
            for i,l in enumerate('ABCD'):
                plt.plot(tout, yout[:,i], label=l)
            plt.legend(loc='best')
            if show:
                plt.show()
            else:
                plt.savefig(__file__[:-2]+'png')


if __name__ == '__main__':
    argh.dispatch_command(main)
