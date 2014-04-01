#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac import ReactionDiffusion, FLAT, CYLINDRICAL, SPHERICAL, Geom_names
from chemreac.integrate import run

def flat_analytic(x, t, D, mu):
    return (4*np.pi*D*t)**-0.5 * np.exp(-(x-mu)**2/(4*D*t))

def spherical_analytic(x, t, D, mu):
    return (4*np.pi*D)**-0.5 * t**-1.5 * np.exp(-(x-mu)**2/(4*D*t))

def cylindrical_analytic(x, t, D, mu):
    return (4*np.pi*D*t)**-1 * np.exp(-(x-mu)**2/(4*D*t))

def main(D=2e-3, t0=1., tend=2., x0=1., xend=2., mu=None, N=1e5, nt=30, geom='f'):
    mu = mu or .5*(x0+xend)
    tout = np.linspace(t0, tend, nt)

    assert geom in 'fcs'
    geom = {'f': FLAT, 'c': CYLINDRICAL, 's': SPHERICAL}[geom]
    print(Geom_names[geom])
    analytic = {
        FLAT: flat_analytic,
        CYLINDRICAL: cylindrical_analytic,
        SPHERICAL: spherical_analytic
    }[geom]

    sys = ReactionDiffusion(
        1, [], [], [], N, D=[D], x=np.linspace(x0, xend, N+1), geom=geom)
    y0 = analytic(sys.x_centers, t0, D, mu)
    t = tout.copy().reshape((nt,1))
    yref = analytic(sys.x_centers, t, D, mu)
    yout, info = run(sys, y0, tout, atol=1e-6, rtol=1e-6)

    # Plot results
    def plot(y, c, ttl=None):
        plt.plot(sys.x_centers, y, c=c)
        plt.xlabel('x / m')
        plt.ylabel('C / M')
        if ttl: plt.title(ttl)

    for i in range(nt):
        c = 1-tout[i]/tend
        c = (1.0-c, .5-c/2, .5-c/2)
        plt.subplot(4,1,1)
        plot(yout[i,:], c, 'Simulation (N={})'.format(sys.N))
        plt.subplot(4,1,2)
        plot(yref[i,:], c, 'Analytic')
        plt.subplot(4,1,3)
        plot((yref[i,:]-yout[i,:])/info['atol'], c,
             'Abs. err. / Abs. tol. (={})'.format(info['atol']))
    plt.subplot(4,1,4)
    plt.plot(tout, np.sum((yref-yout)**2/N, axis=1)**0.5/info['atol'])
    plt.xlabel('Time / s')
    plt.ylabel(r'$\sqrt{\langle E^2 \rangle} / atol$')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
