#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argh
import numpy as np

from chemreac.serialization import load
from chemreac.cpp_chem_wrapper import \
    PyReactionDiffusion as ReactionDiffusion

from chemreac.integrate import run


"""
Demo of chemical reaction diffusion system.
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0

def main(tend=10.0, N=1, nt=50, plot=True):
    sys = load('four_species.json', N=N)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4]+\
                  [0.1, 0.1, 0.1, 0.1]*(sys.N-1))


    t0 = 0.0
    tout, yout, info = run(sys, y0, t0, tend, nt)
    if plot:
        plt.plot(tout, yout[:,0], label='A')
        plt.plot(tout, yout[:,1], label='B')
        plt.plot(tout, yout[:,2], label='C')
        plt.plot(tout, yout[:,3], label='D')

        plt.show()

    # Used for testing
    assert sys.N == 1 # Leave diffusion of of the picture for now
    k1, k2 = sys.k
    A, B, C, D = y0
    ref_f = np.array([-k1*A, k1*A, -2*k2*C*C*B, k2*C*C*B])

    ref_J = np.array([[-k1,         0,           0, 0],
                      [ k1,         0,           0, 0],
                      [  0, -2*k2*C*C, -2*2*k2*C*B, 0],
                      [  0,    k2*C*C,    2*k2*C*B, 0]])

    print(sys, y0, t0, ref_f, ref_J)

if __name__ == '__main__':
    argh.dispatch_command(main)
