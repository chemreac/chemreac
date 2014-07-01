#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import
from future.builtins import *

import argh
import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import run
from chemreac.util.analysis import solver_linear_error

analytic = {
    0: lambda y0, k, t: y0[0] * np.exp(-k[0]*t),
    1: lambda y0, k, t: (y0[1] * np.exp(-k[1] * t) +
                        y0[0] * k[0] / (k[1] - k[0]) *
                        (np.exp(-k[0]*t) -
                         np.exp( - k[1] * t))),
    2: lambda y0, k, t: (y0[2] * np.exp(-k[2] * t) +
                        y0[1] * k[1] / (k[2] - k[1]) *
                        (np.exp(-k[1]*t) -
                         np.exp(-k[2]*t)) +
                        k[1] * k[0] * y0[0] / (k[1] - k[0]) *
                        (1 / (k[2] - k[0]) *
                         (np.exp( - k[0] * t) -
                          np.exp( - k[2] * t)) -
                         1 / (k[2] - k[1]) *
                         (np.exp( - k[1] * t) -
                          np.exp( - k[2] * t))))
}

def get_yref(k, y0, tout):
    coeffs = k + [0]*(3-len(k))
    return np.column_stack([
        analytic[i](y0, coeffs, tout) for i in range(min(3,len(k)+1))
    ])


def integrate_rd(t0=1e-10, tend=2.0, A0=42.42, nt=67, small=20,
                 rates='3.0,4.0', logy=False, logt=False,
                 plot=False, atol='1e-7,1e-6,1e-5', rtol='1e-6',
                 num_jac=False):
    """
    Simplest possible reaction system: 1st order decay:
    A -> B

    Analytic solution through Bateman equation => ensure |k_i - k_j| >> eps
    """
    k = list(map(float, rates.split(',')))
    n = len(k)+1
    if n > 4:
        raise ValueError("Maximum 3 consequtive decays supported at the moment.")

    atol = list(map(float, atol.split(',')))
    if len(atol) == 1:
        atol = atol[0]

    rtol = list(map(float, rtol.split(',')))
    if len(rtol) == 1:
        rtol = rtol[0]

    rd = ReactionDiffusion(
        n, [[i] for i in range(n-1)], [[i] for i in range(1, n)],
        k, logy=logy, logt=logt)

    y0 = np.array([A0] + [10**-small]*(n-1))
    tout = np.linspace(t0, tend, nt)
    yref = get_yref(k, y0, tout)

    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(rd, y, t, atol=atol, rtol=rtol,
                     with_jacobian=not num_jac)
    linC = np.exp(yout) if logy else yout
    linC = linC[:, 0, :]

    if plot:
        nshow=min(n,3)
        try:
            min_atol = min(info['atol'])
        except:
            min_atol = info['atol']

        import matplotlib.pyplot as plt
        c = 'rgb'
        for i, l in enumerate('ABC'[:nshow]):
            plt.subplot(nshow+1, 1, 1)
            plt.plot(tout, linC[:, i], label=l, color=c[i])

            plt.subplot(nshow+1, 1, 2+i)
            try:
                atol = info['atol'][i]
            except:
                atol = info['atol']

            try:
                rtol = info['rtol'][i]
            except:
                rtol = info['rtol']
            plt.plot(tout, (linC[:, i]-yref[:, i])/min_atol, label=l, color=c[i])
            le_l, le_u = solver_linear_error(yout[:, 0, i], atol, rtol,
                                             rd.logy) - linC[:, i]
            plt.fill_between(tout, le_l/min_atol, le_u/min_atol, color=c[i], alpha=0.2)

        plt.subplot(nshow+1, 1, 1)
        plt.title('Concentration vs. time')
        plt.legend(loc='best', prop={'size': 11})
        plt.xlabel('t')
        plt.ylabel('[X]')
        for i in range(nshow):
            plt.subplot(nshow+1, 1, 2+i)
            plt.title('Absolute error in [{}](t) / min(atol)'.format('ABC'[i]))
            plt.legend(loc='best')
            plt.xlabel('t')
            plt.ylabel('|E[{0}]| / {1:7.0g}'.format('ABC'[i], min_atol))
        plt.tight_layout()
        plt.show()

    return linC, yref, rd, info

if __name__ == '__main__':
    argh.dispatch_command(integrate_rd)
