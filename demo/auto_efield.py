#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

from math import log, erf, exp

import argh
import numpy as np

from chemreac import (
    ReactionDiffusion, FLAT, CYLINDRICAL, SPHERICAL, Geom_names
)
from chemreac.integrate import run
from chemreac.util.transforms import sigmoid_algebraic_4

sq2 = 2**0.5
pi = np.pi
sqpi = pi**0.5

def gaussian(x, mu, sigma, logy, logx, geom):
    # Mathematica code
    # $Assumptions = {(sigma | mu) \[Element] Reals, sigma > 0}
    # 1/Integrate[E^(-1/2*((x - mu)/sigma)^2), {x, -Infinity, Infinity}]
    # 1/Integrate[2*pi*x*E^(-1/2*((x - mu)/sigma)^2), {x, 0, Infinity}]
    # 1/Integrate[4*pi*x^2*E^(-1/2*((x - mu)/sigma)^2), {x, 0, Infinity}]
    if geom == FLAT:
        a = 1/sigma/(2*np.pi)**0.5
    elif geom == CYLINDRICAL:
        a = 1/pi/sigma/(2*exp(-mu**2/2/sigma**2)*sigma +\
                        mu*sq2*sqpi*(1 + erf(mu/(sq2*sigma))))
    elif geom == SPHERICAL:
        a = 1/2/pi/sigma/(2*exp(-mu**2/2/sigma**2)*mu*sigma +\
                          sq2*sqpi*(mu**2 + sigma**2)*(1 + erf(mu/sq2/sigma)))
    else:
        raise RuntimeError()
    b = -0.5*((x-mu)/sigma)**2
    if logy:
        return log(a) + b
    else:
        return a*np.exp(b)

def pair_of_gaussians(x, offsets, sigma, logy, logx, geom):
    x = np.exp(x) if logx else x
    xspan = (x[-1] - x[0])
    xl = x[0] + (0.5 + offsets[0])*xspan  # lower
    xu = x[0] + (0.5 + offsets[1])*xspan  # upper
    return (
        gaussian(x, xl, sigma, logy, logx, geom),
        gaussian(x, xu, sigma, logy, logx, geom)
    )


def integrate_rd(D=0., t0=1e-6, tend=7., x0=0.1, xend=1.0, N=256,
                 offset=0.25, mobility=3e-8, nt=25, geom='f',
                 logt=False, logy=False, logx=False, random=False,
                 nstencil=3, lrefl=False, rrefl=False,
                 num_jacobian=False, method='bdf', plot=False,
                 atol=1e-6, rtol=1e-6, random_seed=42, surf_chg=0.0,
                 sigma_q=101):
    if random_seed:
        np.random.seed(random_seed)
    n = 2
    tout = np.linspace(t0, tend, nt)
    geom = {'f': FLAT, 'c': CYLINDRICAL, 's': SPHERICAL}[geom]

    # Setup the grid
    _x0 = log(x0) if logx else x0
    _xend = log(xend) if logx else xend
    x = np.linspace(_x0, _xend, N+1)
    if random:
        x += (np.random.random(N+1)-0.5)*(_xend-_x0)/(N+2)

    # Setup the system
    stoich_reac = []
    stoich_prod = []
    k = []
    bin_k_factor = [[] for _ in range(N)]
    bin_k_factor_span = []

    sys = ReactionDiffusion(
        n, stoich_reac, stoich_prod, k, N,
        D=[D, D],
        z_chg=[1, -1],
        mobility=[mobility, -mobility],
        x=x,
        bin_k_factor=bin_k_factor,
        bin_k_factor_span=bin_k_factor_span,
        geom=geom,
        logy=logy,
        logt=logt,
        logx=logx,
        nstencil=nstencil,
        lrefl=lrefl,
        rrefl=rrefl,
        auto_efield=True,
        surf_chg=surf_chg,
        eps=80.10,  # water at 20 deg C
    )

    # Initial conditions
    s = (xend-x0)/sigma_q
    y0 = np.vstack(pair_of_gaussians(sys.xcenters, [offset, -offset],
                                     s, logy, logx, geom)).transpose()
    if logy:
        y0 = sigmoid_algebraic_4(y0)
    print(y0)
    if plot:
        # Plot initial E-field
        import matplotlib.pyplot as plt
        sys.calc_efield((np.exp(y0) if logy else y0).flatten())
        plt.subplot(4, 1, 3)
        plt.plot(sys.xcenters, sys.efield)
        plt.plot(sys.xcenters, sys.xcenters*0)



    # Run the integration
    t = np.log(tout) if logt else tout
    yout, info = run(sys, y0.flatten(), t,
                     atol=atol, rtol=rtol,
                     with_jacobian=(not num_jacobian), method=method)
    yout = np.exp(yout) if logy else yout

    # Plot results
    if plot:
        def _plot(y, ttl=None,  **kwargs):
            plt.plot(sys.xcenters, y, **kwargs)
            if N < 100:
                plt.vlines(sys.x, 0, np.ones_like(sys.x)*max(y),
                           linewidth=.1, colors='gray')
            plt.xlabel('x / m')
            plt.ylabel('C / M')
            if ttl:
                plt.title(ttl)

        for i in range(nt):
            plt.subplot(4, 1, 1)
            c = 1-tout[i]/tend
            c = (1.0-c, .5-c/2, .5-c/2)
            _plot(yout[i, :, 0], 'Simulation (N={})'.format(sys.N),
                  c=c, label='$z_A=1$' if i==0 else None)
            _plot(yout[i, :, 1], c=c[::-1],
                  label='$z_B=-1$' if i==0 else None)
            plt.legend()

            plt.subplot(4, 1, 2)
            delta_y = yout[i, :, 0] - yout[i, :, 1]
            _plot(delta_y, 'Diff'.format(sys.N),
                  c=[c[2], c[0], c[1]],
                  label='A-B (positive excess)' if i==0 else None)
            plt.legend()
            plt.xlabel('Time / s')
            plt.ylabel(r'C')
        plt.subplot(4, 1, 3)
        plt.plot(sys.xcenters, sys.efield)

        for i in range(3):
            plt.subplot(4, 1, i+1)
            ax = plt.gca()
            for d in (-1, 1):
                plt.plot([x0+(0.5+d*offset)*(xend-x0)]*2,
                         ax.get_ylim(), '--k')
        plt.subplot(4, 1, 4)
        for i in range(n):
            plt.plot(tout, [sys.integrated_conc(yout[j,:,i]) for j in range(nt)])
        plt.tight_layout()
        plt.show()
    return tout, yout, info, sys

if __name__ == '__main__':
    argh.dispatch_command(integrate_rd)
