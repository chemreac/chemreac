#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
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

 $ python analytic_diffusion.py --plot --efield --center 0.5\
 --nstencil 5 --nspecies 3 --geom f

 $ python analytic_diffusion.py --x0 0 --xend 1000 --N 1000 --center 500\
 -D 400 --nstencil 3

Note -D 475

::

 $ python analytic_diffusion.py --x0 0 --xend 1000 --N 1000 --center 500\
 -D 475 --nstencil 7

Still problematic (should not need to be):

::

 $ python analytic_diffusion.py --plot --nstencil 5 --logy --D 0.0005

Here is an example generated by:

::

 $ python analytic_diffusion.py --plot --nstencil 3 --nspecies 2 --geom f\
 --savefig analytic_diffusion.png


.. image:: ../_generated/analytic_diffusion.png


"""

from __future__ import absolute_import, division, print_function


from math import log
from itertools import product

import argh
import numpy as np

from chemreac import ReactionDiffusion
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


def integrate_rd(N=64, geom='f', nspecies=1, nstencil=3,
                 D=2e-3, t0=3.0, tend=7., x0=0.0, xend=1.0, center=None,
                 nt=42, logt=False, logy=False, logx=False,
                 random=False, p=0, a=0.2,
                 linterpol=False, rinterpol=False, ilu_limit=5.0,
                 n_jac_diags=-1, num_jacobian=False,
                 method='bdf', solver='sundials', iter_type='default',
                 linear_solver='default',
                 atol=1e-8, rtol=1e-10,
                 efield=False, random_seed=42, mobility=0.01,
                 plot=False, savefig='None', verbose=False, yscale='linear',
                 vline_limit=100,
                 ):  # remember: anayltic_N_scaling.main kwargs
    if t0 == 0.0:
        raise ValueError("t0==0 => Dirac delta function C0 profile.")
    if random_seed:
        np.random.seed(random_seed)
    # decay = (nspecies > 1)
    # n = 2 if decay else 1
    center = float(center or x0)
    tout = np.linspace(t0, tend, nt)

    assert geom in 'fcs'
    analytic = {
        'f': flat_analytic,
        'c': cylindrical_analytic,
        's': spherical_analytic
    }[geom]

    # Setup the grid
    _x0 = log(x0) if logx else x0
    _xend = log(xend) if logx else xend
    x = np.linspace(_x0, _xend, N+1)
    if random:
        x += (np.random.random(N+1)-0.5)*(_xend-_x0)/(N+2)

    def _k(si):
        return (si+p)*log(a+1)
    k = [_k(i+1) for i in range(nspecies-1)]
    rd = ReactionDiffusion(
        nspecies,
        [[i] for i in range(nspecies-1)],
        [[i+1] for i in range(nspecies-1)],
        k,
        N,
        D=[D]*nspecies,
        z_chg=[1]*nspecies,
        mobility=[mobility]*nspecies,
        x=x,
        geom=geom,
        logy=logy,
        logt=logt,
        logx=logx,
        nstencil=nstencil,
        lrefl=not linterpol,
        rrefl=not rinterpol,
        ilu_limit=ilu_limit,
        n_jac_diags=n_jac_diags
    )

    if efield:
        if geom != 'f':
            raise ValueError("Only analytic sol. for flat drift implemented.")
        rd.efield = _efield_cb(rd.xcenters)

    # Calc initial conditions / analytic reference values
    t = tout.copy().reshape((nt, 1))
    yref = analytic(rd.xcenters, t, D, center, x0, xend,
                    -mobility if efield else 0, logy, logx).reshape(nt, N, 1)

    if nspecies > 1:
        from batemaneq import bateman_parent
        bateman_out = np.array(bateman_parent(k, tout)).T
        terminal = (1 - np.sum(bateman_out, axis=1)).reshape((nt, 1))
        bateman_out = np.concatenate((bateman_out, terminal), axis=1).reshape(
            (nt, 1, nspecies))
        if logy:
            yref = yref + np.log(bateman_out)
        else:
            yref = yref * bateman_out

    # Run the integration
    integr = run(rd, yref[0, ...], tout, atol=atol, rtol=rtol,
                 with_jacobian=(not num_jacobian), method=method,
                 iter_type=iter_type, linear_solver=linear_solver,
                 C0_is_log=logy, solver=solver)
    info = integr.info

    if logy:
        def lin_err(i, j):
            linref = np.exp(yref[i, :, j])
            linerr = np.exp(integr.yout[i, :, j])-linref
            linatol = np.average(yref[i, :, j])
            linrtol = linatol
            return linerr/(linrtol*np.abs(linref)+linatol)

    if logy:
        rmsd = np.sum(lin_err(slice(None), slice(None))**2 / N, axis=1)**0.5
    else:
        rmsd = np.sum((yref-integr.yout)**2 / N, axis=1)**0.5
    ave_rmsd_over_atol = np.average(rmsd, axis=0)/atol

    if verbose:
        # Print statistics
        from pprint import pprint
        pprint(info)
        pprint(ave_rmsd_over_atol)

    # Plot results
    if plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(6, 10))

        # colors: (0.5, 0.5, 0.5), (0.5, 0.5, 1), ...
        base_colors = list(product([.5, 1], repeat=3))[1:-1]

        def color(ci, ti):
            return np.array(base_colors[ci % len(base_colors)])*tout[ti]/tend

        for ti in range(nt):
            plt.subplot(4, 1, 1)
            for si in range(nspecies):
                plt.plot(rd.xcenters, integr.Cout[ti, :, si], c=color(si, ti),
                         label=None if ti < nt - 1 else rd.substance_names[si])

            plt.subplot(4, 1, 2)
            for si in range(nspecies):
                plt.plot(rd.xcenters, np.exp(yref[ti, :, si]) if logy
                         else yref[ti, :, si], c=color(si, ti))

            plt.subplot(4, 1, 3)
            if logy:
                for si in range(nspecies):
                    plt.plot(rd.xcenters, lin_err(ti, si)/atol,
                             c=color(si, ti))
            else:
                for si in range(nspecies):
                    plt.plot(
                        rd.xcenters,
                        (yref[ti, :, si] - integr.yout[ti, :, si])/atol,
                        c=color(si, ti))

        if N < vline_limit:
            for idx in range(1, 4):
                plt.subplot(4, 1, idx)
                for bi in range(N):
                    plt.axvline(rd.x[bi], color='gray')

        plt.subplot(4, 1, 1)
        plt.title('Simulation (N={})'.format(rd.N))
        plt.xlabel('x / m')
        plt.ylabel('C / M')
        plt.gca().set_yscale(yscale)
        plt.legend()

        plt.subplot(4, 1, 2)
        plt.title('Analytic solution')
        plt.gca().set_yscale(yscale)

        plt.subplot(4, 1, 3)
        plt.title('Linear rel. error / Abs. tol. (={})'.format(atol))

        plt.subplot(4, 1, 4)
        plt.title('RMS error vs. time'.format(atol))
        tspan = [tout[0], tout[-1]]
        for si in range(nspecies):
            plt.plot(tout, rmsd[:, si] / atol, c=color(si, -1))
            plt.plot(tspan, [ave_rmsd_over_atol[si]]*2,
                     c=color(si, -1), ls='--')

        plt.xlabel('Time / s')
        plt.ylabel(r'$\sqrt{\langle E^2 \rangle} / atol$')
        plt.tight_layout()
        save_and_or_show_plot(savefig=savefig)

    return tout, integr.yout, info, ave_rmsd_over_atol, rd, rmsd


if __name__ == '__main__':
    argh.dispatch_command(integrate_rd, output_file=None)
