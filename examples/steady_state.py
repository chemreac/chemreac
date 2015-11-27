#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function
from future.builtins import *

from math import log

import argh
import numpy as np

from chemreac import (
    ReactionDiffusion, FLAT, CYLINDRICAL, SPHERICAL
)
from chemreac.integrate import run
from chemreac.util.plotting import plot_solver_linear_error


def efield_cb(x, logx=False):
    """
    Returns a flat efield (-1)
    """
    return -np.ones_like(x)


def y0_flat_cb(x, logx=False):
    xc = x[:-1] + np.diff(x)/2
    if logx:
        x, xc = map(np.exp, (x, xc))
    return 17 - 11*(xc-x[0])/(x[-1]-x[0])


def y0_cylindrical_cb(x, logx=False):
    xc = x[:-1] + np.diff(x)/2
    if logx:
        x, xc = map(np.exp, (x, xc))
    return 17 - np.log((xc-x[0])/(x[-1]-x[0]))


def y0_spherical_cb(x, logx=False):
    xc = x[:-1] + np.diff(x)/2
    if logx:
        x, xc = map(np.exp, (x, xc))
    return 3 + 0.1/((xc-x[0])/(x[-1]-x[0]))


def integrate_rd(D=2e-3, t0=3., tend=7., x0=0.0, xend=1.0, mu=None, N=32,
                 nt=25, geom='f', logt=False, logy=False, logx=False,
                 random=False, nstencil=3, lrefl=False, rrefl=False,
                 num_jacobian=False, method='bdf', plot=False,
                 atol=1e-6, rtol=1e-6, efield=False, random_seed=42,
                 verbose=False):
    if random_seed:
        np.random.seed(random_seed)
    n = 1
    mu = float(mu or x0)
    tout = np.linspace(t0, tend, nt)

    assert geom in 'fcs'
    geom = {'f': FLAT, 'c': CYLINDRICAL, 's': SPHERICAL}[geom]

    # Setup the grid
    _x0 = log(x0) if logx else x0
    _xend = log(xend) if logx else xend
    x = np.linspace(_x0, _xend, N+1)
    if random:
        x += (np.random.random(N+1)-0.5)*(_xend-_x0)/(N+2)

    mob = 0.3
    # Initial conditions
    y0 = {
        FLAT: y0_flat_cb,
        CYLINDRICAL: y0_cylindrical_cb,
        SPHERICAL: y0_spherical_cb
    }[geom](x, logx)

    # Setup the system
    stoich_active = []
    stoich_prod = []
    k = []

    assert not lrefl
    assert not rrefl
    # bin_k_factor = [[] for _ in range(N)]
    # bin_k_factor_span = []
    # if lrefl:
    #     # isolated sys is not stationary for linear conc profile with finite
    #     # slope, hence we add production reaction at left boundary
    #     assert nstencil == 3  # k is derived for parabola through -x0, x0, x1
    #     stoich_active.append([0])
    #     stoich_prod.append([0, 0])
    #     for i in range(N):
    #         bin_k_factor[i].append(1 if i == 0 else 0)
    #     bin_k_factor_span.append(1)
    #     x0 = (x[0]+x[1])/2 - x[0]
    #     x1 = (x[1]+x[2])/2 - x[0]
    #     C = {FLAT: 2, CYLINDRICAL: 4, SPHERICAL: 6}[geom]
    #     k.append(C*D*(y0[0]-y0[1])/(y0[0]*(x1**2 - x0**2)))
    # if rrefl:
    #     # for same reason as above, a consumption reaction is added at right
    #     # boundary
    #     assert nstencil == 3
    #     stoich_active.append([0])
    #     stoich_prod.append([])
    #     for i in range(N):
    #         bin_k_factor[i].append(1 if i == (N-1) else 0)
    #     bin_k_factor_span.append(1)
    #     xNm1 = (x[N-1]+x[N])/2-x[N]
    #     xNm2 = (x[N-2]+x[N-1])/2-x[N]
    #     C = {FLAT: 2, CYLINDRICAL: 4, SPHERICAL: 6}[geom]
    #     k.append(C*D*(y0[N-2]-y0[N-1])/(y0[N-1]*(xNm2**2 - xNm1**2)))

    rd = ReactionDiffusion(
        n, stoich_active, stoich_prod, k, N,
        D=[D],
        z_chg=[1],
        mobility=[mob],
        x=x,
        # bin_k_factor=bin_k_factor,
        # bin_k_factor_span=bin_k_factor_span,
        geom=geom,
        logy=logy,
        logt=logt,
        logx=logx,
        nstencil=nstencil,
        lrefl=lrefl,
        rrefl=rrefl,
    )

    if efield:
        if geom != FLAT:
            raise ValueError("Only analytic sol. for flat drift implemented.")
        rd.efield = efield_cb(rd.xcenters, logx)

    # Analytic reference values
    t = tout.copy().reshape((nt, 1))
    Cref = np.repeat(y0[np.newaxis, :, np.newaxis], nt, axis=0)
    if efield:
        Cref += t.reshape((nt, 1, 1))*mob

    # Run the integration
    integr = run(rd, y0, tout, atol=atol, rtol=rtol,
                 with_jacobian=(not num_jacobian), method=method)
    Cout, info = integr.Cout, integr.info

    if verbose:
        print(info)

    def lin_err(i=slice(None), j=slice(None)):
        return integr.Cout[i, :, j] - Cref[i, :, j]

    rmsd = np.sum(lin_err()**2 / N, axis=1)**0.5
    ave_rmsd_over_atol = np.average(rmsd, axis=0)/info['atol']

    # Plot results
    if plot:
        import matplotlib.pyplot as plt

        def _plot(y, c, ttl=None, apply_exp_on_y=False):
            plt.plot(rd.xcenters, np.exp(y) if apply_exp_on_y else y, c=c)
            if N < 100:
                plt.vlines(rd.x, 0, np.ones_like(rd.x)*max(y), linewidth=.1,
                           colors='gray')
            plt.xlabel('x / m')
            plt.ylabel('C / M')
            if ttl:
                plt.title(ttl)

        for i in range(nt):
            c = 1-tout[i]/tend
            c = (1.0-c, .5-c/2, .5-c/2)  # over time: dark red -> light red

            plt.subplot(4, 1, 1)
            _plot(Cout[i, :, 0], c, 'Simulation (N={})'.format(rd.N),
                  apply_exp_on_y=logy)

            plt.subplot(4, 1, 2)
            _plot(Cref[i, :, 0], c, 'Analytic', apply_exp_on_y=logy)

            ax_err = plt.subplot(4, 1, 3)
            plot_solver_linear_error(integr, Cref, ax_err, ti=i,
                                     bi=slice(None),
                                     color=c, fill=(i == 0))
            plt.title('Linear rel error / Log abs. tol. (={})'.format(
                      info['atol']))

        plt.subplot(4, 1, 4)
        tspan = [tout[0], tout[-1]]
        plt.plot(tout, rmsd[:, 0] / info['atol'], 'r')
        plt.plot(tspan, [ave_rmsd_over_atol[0]]*2, 'r--')

        plt.xlabel('Time / s')
        plt.ylabel(r'$\sqrt{\langle E^2 \rangle} / atol$')
        plt.tight_layout()
        plt.show()
    return tout, Cout, info, ave_rmsd_over_atol, rd


if __name__ == '__main__':
    argh.dispatch_command(integrate_rd, output_file=None)
