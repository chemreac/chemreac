#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

"""
A     ->  B      k1
    B ->  C      k2
A + B ->  B + C  k3

dA/dt = -k1*A          - k3*A*B
dB/dt =  k1*A - k2*B
dC/dt =         k2*B   + k3*A*B

A+B+C = A0 + B0 + C0 => sum of derivatives = 0 (already satisfied)

Steady state assumption for B (A0 >> B0):
B = k1*A/k2
dA/dt = -k1*(A + k3/k2*A**2)
log(A) - log(k3/k2*A + 1) = -k1*t + (log(A0) - log(k3/k2*A0 + 1))
A/(k3/k2*A + 1) = f(t)
{f(t) = A0/(k3/k2*A0 + 1)*exp(-k1*t)}
A(1 - f(t)*k3/k2) = f(t)
A = f(t)/(1 - f(t)*k3/k2) [*]
A = 1/(1/f(t) - k3/k2)
A = 1/((k3/k2*A0 + 1)/A0*exp(k1*t) - k3/k2)

[*] if k3/k2 << 1:
A = f(t)
"""

import argh
import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import run
from chemreac.util.plotting import plot_C_vs_t_in_bin


def integrate_rd(tend=1.0, k1=7e-1, k2=3e2, k3=7.0,
                 A0=1.0, B0=0.0, C0=0.0, plot=False, nt=1024):
    """
    Runs integration and (optionally) generates plots.
    """
    def f(t):
        return A0/(k3/k2*A0 + 1)*np.exp(-k1*t)

    y = [A0, B0, C0]
    k = [k1, k2, k3]
    rd = ReactionDiffusion(
        3, [[0], [1], [0, 1]], [[1], [2], [1, 2]], k)
    t = np.linspace(0, tend, nt)
    yout, info = run(rd, y, t)
    A_ssB = 1/(1/f(t) - k3/k2)
    A_ssB_2fast = f(t)

    if plot:
        import matplotlib.pyplot as plt
        ax = plt.subplot(3, 1, 1)
        plot_C_vs_t_in_bin(rd, t, yout, ax=ax)
        plt.subplot(3, 1, 2)
        plt.plot(t, yout[:, 0, 0] - A_ssB,
                 label="Abs. err. in A assuming steady state of B")
        plt.legend(loc='best')
        plt.subplot(3, 1, 3)
        plt.plot(t, yout[:, 0, 0] - A_ssB_2fast,
                 label="Abs. err. in A when also assuming k2 >> k3")
        plt.legend(loc='best')
        plt.show()
    ydot = lambda x, y: (-k1*y[0] - k3*y[0]*y[1], k1*y[0] - k2*y[1],
                         k2*y[1] + k3*y[0]*y[1])
    return t, yout, A_ssB, A_ssB_2fast, ydot


if __name__ == '__main__':
    argh.dispatch_command(integrate_rd)
