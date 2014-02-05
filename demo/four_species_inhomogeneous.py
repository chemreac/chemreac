#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

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

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0


def main(tend=10.0, N=10, nt=50, plot=True, spy=False, mode=None):
    mod1 = lambda x: x/(x**2+1)
    sys = load('four_species.json', N=N, x=N,
               bin_k_factor=[[mod1(x)] for x in range(N)], # decay A->B is modulated with x
               bin_k_factor_span=[1]
    )

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    y0 = np.concatenate([y0/(i+1)*(0.25*i**2+1) for i in range(N)])

    t0 = 0.0
    tout, yout, info = run(sys, y0, t0, tend, nt)
    if plot:
        for i,l in enumerate('ABCD'):
            plt.plot(tout, yout[:,i], label=l)
        plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
