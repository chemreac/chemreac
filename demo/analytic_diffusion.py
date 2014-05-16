#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
 $ python analytic_diffusion.py --x0 0 --xend 1000 --N 1000 --mu 500 -D 400 \
     --nstencil 3

Note -D 475

 $ python analytic_diffusion.py --x0 0 --xend 1000 --N 1000 --mu 500 -D 475 \
     --nstencil 7
"""

from __future__ import print_function, division, absolute_import

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac import (
    ReactionDiffusion, FLAT, CYLINDRICAL, SPHERICAL, Geom_names
)
from chemreac.integrate import run


def flat_analytic(x, t, D, mu):
    return (4*np.pi*D*t)**-0.5 * np.exp(-(x-mu)**2/(4*D*t))


def spherical_analytic(x, t, D, mu):
    return (4*np.pi*D)**-0.5 * t**-1.5 * np.exp(-(x-mu)**2/(4*D*t))


def cylindrical_analytic(x, t, D, mu):
    return (4*np.pi*D*t)**-1 * np.exp(-(x-mu)**2/(4*D*t))


def main(D=2e-3, t0=3., tend=7., x0=0.0, xend=1.0, mu=None, N=2048, nt=30,
         geom='f', logt=False, logy=False, random=False, k=0.0, nstencil=3,
         linterpol=False, rinterpol=False, num_jacobian=False, method='bdf',
         scale_x=False):
    decay = (k != 0.0)
    n = 2 if decay else 1
    mu = float(mu or x0)
    tout = np.linspace(t0, tend, nt)

    assert geom in 'fcs'
    geom = {'f': FLAT, 'c': CYLINDRICAL, 's': SPHERICAL}[geom]
    print(Geom_names[geom])
    analytic = {
        FLAT: flat_analytic,
        CYLINDRICAL: cylindrical_analytic,
        SPHERICAL: spherical_analytic
    }[geom]

    # Setup the system
    # x = np.logspace(np.log(x0), np.log(xend), N+1, base=np.exp(1))
    x = np.linspace(x0, xend, N+1)

    if random:
        x += (np.random.random(N+1)-0.5)*(xend-x0)/(N+2)
    sys = ReactionDiffusion(
        2 if decay else 1,
        [[0]] if decay else [],
        [[1]] if decay else [],
        [k] if decay else [],
        N,
        D=[D]*2 if decay else [D],
        x=x,
        geom=geom,
        logy=logy,
        logt=logt,
        nstencil=nstencil,
        lrefl=not linterpol,
        rrefl=not rinterpol,
        xscale=1/(x[1]-x[0]) if scale_x else 1.0
    )

    # Calc initial conditions / analytic reference values
    t = tout.copy().reshape((nt, 1))
    yref = (xend-x0)*analytic(sys.xcenters, t, D, mu).reshape(nt, N, 1)
    if decay:
        yref = np.concatenate((yref, yref), axis=2)
        yref[:, :, 0] *= np.exp(-k*t)
        yref[:, :, 1] *= 1-np.exp(-k*t)

    y0 = yref[0, ...].flatten()

    # Run the integration
    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(sys, y, t, atol=1e-8, rtol=1e-6,
                     with_jacobian=(not num_jacobian), method=method)
    yout = np.exp(yout) if logy else yout
    print(info)

    # Plot results
    def plot(y, c, ttl=None):
        plt.plot(sys.xcenters, y, c=c)
        if N < 100:
            plt.vlines(sys.x, 0, np.ones_like(sys.x)*y0[0], linewidth=.1,
                       colors='gray')
        plt.xlabel('x / m')
        plt.ylabel('C / M')
        if ttl:
            plt.title(ttl)

    for i in range(nt):
        c = 1-tout[i]/tend
        c = (1.0-c, .5-c/2, .5-c/2)

        plt.subplot(4, 1, 1)
        plot(yout[i, :, 0], c, 'Simulation (N={})'.format(sys.N))
        if decay:
            plot(yout[i, :, 1], c[::-1])

        plt.subplot(4, 1, 2)
        plot(yref[i, :, 0], c, 'Analytic')
        if decay:
            plot(yref[i, :, 1], c[::-1])

        plt.subplot(4, 1, 3)
        plot((yref[i, :, 0]-yout[i, :, 0])/info['atol'], c,
             'Abs. err. / Abs. tol. (={})'.format(info['atol']))
        if decay:
            plot((yref[i, :, 1]-yout[i, :, 1])/info['atol'], c[::-1])

    plt.subplot(4, 1, 4)
    plt.plot(tout, np.sum((yref[:, :, 0] - yout[:, :, 0])**2 / N,
                          axis=1)**0.5 / info['atol'], 'r')
    if decay:
        plt.plot(tout, np.sum((yref[:, :, 1]-yout[:, :, 1])**2/N,
                              axis=1)**0.5/info['atol'], 'b')
    plt.xlabel('Time / s')
    plt.ylabel(r'$\sqrt{\langle E^2 \rangle} / atol$')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
