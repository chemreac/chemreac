#!/usr/bin/env python
# -*- coding: utf-8 -*-

u"""
Analytic diffusion
------------------

:download:`examples/analytic_diffusion.py` models a diffusion process
and reports the error from the model integration by comparison to the
analytic solution (intial concentrations are taken from Green's
function expressions for respective geometry).

::

 $ python analytic_diffusion.py --help

.. exec::
   echo "::\\n\\n"
   python examples/examples/analytic_diffusion.py --help | sed "s/^/   /"

::

 $ python analytic_diffusion.py --plot --efield --mu 0.5 --nstencil 5 --k 0.1\
 --geom f

 $ python analytic_diffusion.py --x0 0 --xend 1000 --N 1000 --mu 500 -D 400\
 --nstencil 3

Note -D 475

::

 $ python analytic_diffusion.py --x0 0 --xend 1000 --N 1000 --mu 500 -D 475\
  --nstencil 7

Still problematic (should not need to be):

::

 $ python analytic_diffusion.py --plot --nstencil 5 --logy --D 0.0005

Here is an example generated by:

::

 $ python analytic_diffusion.py --plot --nstencil 3 --k 0.1 --geom f\
 --savefig analytic_diffusion.png


.. image:: ../_generated/analytic_diffusion.png


"""

from __future__ import absolute_import, division, print_function
from future.builtins import *


from math import log

import argh
import numpy as np

from chemreac import (
    ReactionDiffusion, FLAT, CYLINDRICAL, SPHERICAL
)
from chemreac.integrate import run
from chemreac.util.plotting import save_and_or_show_plot


def flat_analytic(x, t, D, mu, x0, xend, v, logy=False, logx=False):
    r"""
    Evaluates the Green's function:

    .. math ::

        c(x, t) = \frac{x_{end}-x_{0}}{\sqrt{4 \pi D t}} \
            e^{-\frac{(x - \mu - vt)^2}{4Dt}}

    which satisfies:

    .. math ::

        \frac{\partial c(x, t)}{\partial t} = D\nabla^2 \
            c(x, t)  - \vec{v} \cdot c(x, t)

    where :math:`\nabla` in cartesian coordinates with planar (yz-plane)
    symmetry is:

    .. math ::

        \nabla = \frac{\partial}{\partial x}

    and where :math:`\nabla^2` :

    .. math ::

        \nabla^2 = \frac{\partial^2}{\partial x^2}


    """
    x = np.exp(x) if logx else x
    a = (4*np.pi*D*t)**-0.5
    b = -(x-mu-v*t)**2/(4*D*t)
    if logy:
        return np.log(a) + b + log(xend-x0)
    else:
        return a*np.exp(b)*(xend-x0)


def cylindrical_analytic(x, t, D, mu, x0, xend, v, logy=False, logx=False):
    r"""
    Evaluates the Green's function:

    .. math ::

        c(x, t) = \frac{x_{end}-x_{0}}{4 \pi D t} \
            e^{-\frac{(x - \mu - vt)^2}{4Dt}}

    which satisfies:

    .. math ::

        \frac{\partial c(x, t)}{\partial t} = D\nabla^2 \
            c(x, t) - \vec{v} \cdot c(x, t)

    where :math:`\nabla` in cylindrical coordinates with axial symmetry is:

    .. math ::

        \nabla = \frac{1}{x}

    and where :math:`\nabla^2` is:

    .. math ::

        \nabla^2 = \frac{1}{x} \frac{\partial}{\partial x}
            \left( x \frac{\partial}{\partial x} \right)


    """
    x = np.exp(x) if logx else x
    a = (4*np.pi*D*t)**-1
    b = -(x-mu-v*t)**2/(4*D*t)
    if logy:
        return np.log(a) + b + log(xend-x0)
    else:
        return a*np.exp(b)*(xend-x0)


def spherical_analytic(x, t, D, mu, x0, xend, v, logy=False, logx=False):
    r"""
    Evaluates the Green's function:

    .. math ::

        c(x, t) = \frac{x_{end}-x_{0}}{\sqrt{4 \pi D t^3}} \
            e^{-\frac{(x - \mu - vt)^2}{4Dt}}

    which satisfies:

    .. math ::

        \frac{\partial c(x, t)}{\partial t} = D\nabla^2
            c(x, t) - \vec{v} \cdot c(x, t)

    where :math:`\nabla` in spherical coordinates for a isotropic system is:

    .. math ::

        \nabla = \frac{\partial}{\partial x}

    and where :math:`\nabla^2` is:

    .. math ::

        \nabla^2 = \frac{1}{x^2} \frac{\partial}{\partial x}
            \left( x^2 \frac{\partial}{\partial x} \right)


    """
    x = np.exp(x) if logx else x
    a = (4*np.pi*D)**-0.5 * t**-1.5
    b = -(x-mu-v*t)**2/(4*D*t)
    if logy:
        return np.log(a) + b + log(xend-x0)
    else:
        return a*np.exp(b)*(xend-x0)


def _efield_cb(x):
    """
    Returns a flat efield (-1)
    """
    return -np.ones_like(x)


def integrate_rd(D=2e-3, t0=3., tend=7., x0=0.0, xend=1.0, mu=None, N=64,
                 nt=42, geom='f', logt=False, logy=False, logx=False,
                 random=False, k=0.0, nstencil=3, linterpol=False,
                 rinterpol=False, num_jacobian=False, method='bdf',
                 plot=False, atol=1e-6, rtol=1e-6, efield=False,
                 random_seed=42, savefig='None', verbose=False):
    if t0 == 0.0:
        raise ValueError("t0==0 => Dirac delta function C0 profile.")
    if random_seed:
        np.random.seed(random_seed)
    decay = (k != 0.0)
    n = 2 if decay else 1
    mu = float(mu or x0)
    tout = np.linspace(t0, tend, nt)

    assert geom in 'fcs'
    geom = {'f': FLAT, 'c': CYLINDRICAL, 's': SPHERICAL}[geom]
    analytic = {
        FLAT: flat_analytic,
        CYLINDRICAL: cylindrical_analytic,
        SPHERICAL: spherical_analytic
    }[geom]

    # Setup the grid
    _x0 = log(x0) if logx else x0
    _xend = log(xend) if logx else xend
    x = np.linspace(_x0, _xend, N+1)
    if random:
        x += (np.random.random(N+1)-0.5)*(_xend-_x0)/(N+2)

    sys = ReactionDiffusion(
        n,
        [[0]] if decay else [],
        [[1]] if decay else [],
        [k] if decay else [],
        N,
        D=[D]*(2 if decay else 1),
        z_chg=[1]*(2 if decay else 1),
        mobility=[0.01]*(2 if decay else 1),
        x=x,
        geom=geom,
        logy=logy,
        logt=logt,
        logx=logx,
        nstencil=nstencil,
        lrefl=not linterpol,
        rrefl=not rinterpol
    )

    if efield:
        if geom != FLAT:
            raise ValueError("Only analytic sol. for flat drift implemented.")
        sys.efield = _efield_cb(sys.xcenters)

    # Calc initial conditions / analytic reference values
    t = tout.copy().reshape((nt, 1))
    yref = analytic(sys.xcenters, t, D, mu, x0, xend,
                    0.01 if efield else 0, logy, logx).reshape(nt, N, 1)

    if decay:
        yref = np.concatenate((yref, yref), axis=2)
        if logy:
            yref[:, :, 0] += -k*t
            yref[:, :, 1] += np.log(1-np.exp(-k*t))
        else:
            yref[:, :, 0] *= np.exp(-k*t)
            yref[:, :, 1] *= 1-np.exp(-k*t)

    # Run the integration
    integr = run(sys, yref[0, ...], tout, atol=atol, rtol=rtol,
                 with_jacobian=(not num_jacobian), method=method,
                 C0_is_log=logy)
    yout, info = integr.yout, integr.info

    if logy:
        def lin_err(i, j):
            linref = np.exp(yref[i, :, j])
            linerr = np.exp(yout[i, :, j])-linref
            linatol = np.average(yref[i, :, j])
            linrtol = linatol
            return linerr/(linrtol*np.abs(linref)+linatol)

    if logy:
        rmsd = np.sum(lin_err(slice(None), slice(None))**2 / N, axis=1)**0.5
    else:
        rmsd = np.sum((yref-yout)**2 / N, axis=1)**0.5
    ave_rmsd_over_atol = np.average(rmsd, axis=0)/info['atol']

    if verbose:
        # Print statistics
        from pprint import pprint
        pprint(info)
        pprint(ave_rmsd_over_atol)

    # Plot results
    if plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 10))

        def _plot(y, c, ttl=None, apply_exp_on_y=False, vlines=False):
            plt.plot(sys.xcenters, np.exp(y) if apply_exp_on_y else y,
                     c=c)
            if vlines:
                plt.vlines(sys.x, 0, np.ones_like(sys.x)*max(y),
                           linewidth=1, colors='gray')
            plt.xlabel('x / m')
            plt.ylabel('C / M')
            if ttl:
                plt.title(ttl)

        for i in range(nt):
            c = 1-tout[i]/tend
            c = (1.0-c, .5-c/2, .5-c/2)

            plt.subplot(4, 1, 1)
            _plot(yout[i, :, 0], c, 'Simulation (N={})'.format(sys.N),
                  apply_exp_on_y=logy, vlines=(i == 0 and N < 100))
            if decay:
                _plot(yout[i, :, 1], c[::-1], apply_exp_on_y=logy)

            plt.subplot(4, 1, 2)
            _plot(yref[i, :, 0], c, 'Analytic', apply_exp_on_y=logy,
                  vlines=(i == 0 and N < 100))
            if decay:
                _plot(yref[i, :, 1], c[::-1], apply_exp_on_y=logy)

            plt.subplot(4, 1, 3)
            if logy:
                _plot(lin_err(i, 0)/info['atol'], c,
                      'Linear rel error / Log abs. tol. (={})'.format(
                          info['atol']), vlines=(i == nt-1 and N < 100))
                if decay:
                    _plot(lin_err(i, 1)/info['atol'], c[::-1])
            else:
                _plot((yref[i, :, 0]-yout[i, :, 0])/info['atol'], c,
                      'Abs. err. / Abs. tol. (={})'.format(info['atol']),
                      vlines=(i == nt-1 and N < 100))
                if decay:
                    _plot((yref[i, :, 1]-yout[i, :, 1])/info['atol'],
                          c[::-1])

        plt.subplot(4, 1, 4)
        tspan = [tout[0], tout[-1]]
        plt.plot(tout, rmsd[:, 0] / info['atol'], 'r')
        plt.plot(tspan, [ave_rmsd_over_atol[0]]*2, 'r--')
        if decay:
            plt.plot(tout, rmsd[:, 1]/info['atol'], 'b')
            plt.plot(tspan, [ave_rmsd_over_atol[1]]*2, 'b--')

        plt.xlabel('Time / s')
        plt.ylabel(r'$\sqrt{\langle E^2 \rangle} / atol$')
        plt.tight_layout()
        save_and_or_show_plot(savefig=savefig)

    return tout, yout, info, ave_rmsd_over_atol, sys


if __name__ == '__main__':
    argh.dispatch_command(integrate_rd, output_file=None)
