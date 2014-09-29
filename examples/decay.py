#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

from itertools import chain

import argh
import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import run
from chemreac.util.analysis import solver_linear_error, suggest_t0

"""
Motivation for sigmoid damped exp(); vary tend: 5, 700, 1700.
Never mind 700 not being correctly represented, the problem
is 1700 completely ruining the integration (NaN's due to overflow).

 $ python decay.py --plot --rates 1.0 --logy --logt --rtol 1e-13 --atol 1e-6 \
      --scale-err 100.0 --plotlogy --nt 1024 --tend 1700
"""

analytic = {
    0: lambda y0, k, t: y0[0] * np.exp(-k[0]*t),
    1: lambda y0, k, t: (y0[1] * np.exp(-k[1] * t) +
                         y0[0] * k[0] / (k[1] - k[0]) *
                         (np.exp(-k[0]*t) -
                          np.exp(-k[1]*t))),
    2: lambda y0, k, t: (y0[2] * np.exp(-k[2] * t) +
                         y0[1] * k[1] / (k[2] - k[1]) *
                         (np.exp(-k[1]*t) -
                          np.exp(-k[2]*t)) +
                         k[1] * k[0] * y0[0] / (k[1] - k[0]) *
                         (1 / (k[2] - k[0]) *
                          (np.exp(-k[0]*t) -
                           np.exp(-k[2]*t)) -
                          1 / (k[2] - k[1]) *
                          (np.exp(-k[1]*t) -
                           np.exp(-k[2]*t))))
}


def get_linCref(k, y0, tout):
    coeffs = k + [0]*(3-len(k))
    return np.column_stack([
        analytic[i](y0, coeffs, tout) for i in range(
            min(3, len(k)+1))])


def integrate_rd(tend=2.0, A0=42.42, nt=67, t0=0.0,
                 rates='3.0,4.0', logy=False, logt=False,
                 small=20, plot=False, atol='1e-7,1e-6,1e-5', rtol='1e-6',
                 num_jac=False, scale_err=1.0, verbose=False,
                 plotlogy=False, plotlogt=False):
    """
    Simplest possible reaction system: 1st order decay:
    A -> B

    Analytic solution through Bateman equation => ensure |k_i - k_j| >> eps
    """
    k = list(map(float, rates.split(',')))
    n = len(k)+1
    if n > 4:
        raise ValueError("Max 3 consequtive decays supported at the moment.")

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
    if t0 == 0.0 and (logt or logy):
        t0_set = True
        t0 = suggest_t0(rd, y0)
        tend += t0  # same total time
    else:
        t0_set = False
    y = np.log(y0) if logy else y0
    tout = np.linspace(t0, tend, nt)
    t = np.log(tout) if logt else tout
    yout, info = run(rd, y, t, atol=atol, rtol=rtol,
                     with_jacobian=not num_jac)
    linC = np.exp(yout) if logy else yout
    linCref = get_linCref(k, y0, tout - tout[0]).reshape((nt, 1, n))

    if plot:
        print(info)
        print('rate: ', k)
        if t0_set:
            print("t0 = {}".format(t0))
        nshow = min(n, 3)
        try:
            min_atol = min(info['atol'])
        except:
            min_atol = info['atol']

        import matplotlib.pyplot as plt
        c = 'rgb'
        for i, l in enumerate('ABC'[:nshow]):
            ax = plt.subplot(nshow+1, 1, 1)
            if plotlogy:
                    ax.set_yscale('log')
            if plotlogt:
                    ax.set_xscale('log')
            ax.plot(tout, linC[:, 0, i], label=l, color=c[i])

            ax = plt.subplot(nshow+1, 1, 2+i)
            if plotlogy:
                    ax.set_yscale('symlog')  # abs error might be < 0
            if plotlogt:
                    ax.set_xscale('log')
            ax.plot(tout, (linC[:, 0, i]-linCref[:, 0, i])/min_atol,
                    label=l, color=c[i])

            try:
                atol = info['atol'][i]
            except:
                atol = info['atol']

            try:
                rtol = info['rtol'][i]
            except:
                rtol = info['rtol']

            le_l, le_u = solver_linear_error(
                yout[:, 0, i], rtol, atol, rd.logy, scale_err=scale_err)
            plt.fill_between(tout, (le_l - linC[:, 0, i])/min_atol,
                             (le_u - linC[:, 0, i])/min_atol,
                             color=c[i], alpha=0.2)

            # Print indices and values of violations of (scaled) error bounds
            def _print(violation):
                print(violation)
                print(le_l[violation],
                      linCref[violation, 0, i],
                      le_u[violation])
            l_viols = np.where(le_l > linCref[:, 0, i])[0]
            u_viols = np.where(le_u < linCref[:, 0, i])[0]
            if len(l_viols) > 0 or len(u_viols) > 0:
                print("Outside error bounds for rtol, atol:", rtol, atol)
            if verbose:
                for violation in chain(l_viols, u_viols):
                    _print(violation)

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

    return yout, linCref, rd, info

if __name__ == '__main__':
    argh.dispatch_command(integrate_rd, output_file=None)
