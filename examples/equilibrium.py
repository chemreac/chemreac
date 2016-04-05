#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Equilibrium
-----------

:download:`examples/equilibrium.py` demonstrates how
scaling can be used together with tolerances to achieve
desired accuracy from the numerical integration.

We will consider the transient towards an equilibrium
for a dimerization:

.. math ::

    A + &B &\\overset{k_f}{\\underset{k_b}{\\rightleftharpoons}}  C

The analytic solution is (its derivation is left as an exercise):

.. math::

    A(t) &= A_0 - x(t) \\\\
    B(t) &= B_0 - x(t) \\\\
    C(t) &= C_0 + x(t) \\\\
    x(t) &= \\frac{(U-b)(U+b)(e^{Ut}-1)}{2k_f(Ue^{Ut} + U - qe^{Ut} +  q)} \\\\

where

.. math::

    U    &= \\sqrt{A^2k_f^2 + 2ABk_f^2 - 2Ak_bk_f + B^2k_f^2 -
                2Bk_bk_f + 4Ck_bk_f + k_b^2} \\\\
    q    &= Ak_f + Bk_f - k_b


::

 $ python equilibrium.py --help

.. exec::
   echo "::\\\\n\\\\n"
   python examples/examples/equilibrium.py --help | sed "s/^/   /"


Here is an example generated by:

::

 $ python equilibrium.py --plot --savefig equilibrium.png


.. image:: ../_generated/equilibrium.png


If concentrations are far from 1 (and below abstol) the accuracy of
the numerical solution will be very poor:

::

 $ python equilibrium.py --A0 1.0 --B0 1e-10 --C0 1e-30 --kf 10 --kb 1 --t0 0\
 --tend 5 --plot --plotlogy --plotlogt --savefig equilibrium_unscaled.png


.. image:: ../_generated/equilibrium_unscaled.png


But by scaling the concentrations so that the smallest is well above the
absolute tolerance we can get accurate results:

::

 $ python equilibrium.py --scaling 1e10 --A0 1.0 --B0 1e-10 --C0 1e-30 --kf 10\
 --kb 1 --t0 0 --tend 5 --plot --plotlogy --plotlogt --savefig\
 equilibrium_scaled.png


.. image:: ../_generated/equilibrium_scaled.png

"""

from __future__ import absolute_import, division, print_function
import argh
import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import Integration
from chemreac.units import (
    SI_base_registry, second, molar, to_unitless, get_derived_unit
)
from chemreac.util.plotting import (
    save_and_or_show_plot, plot_solver_linear_error,
    plot_solver_linear_excess_error
)


def _algebraic_sigmoid(x, power, limit):
    # Avoid overflow in exp()
    return x/((x/limit)**power+1)**(-1/power)


def analytic_x(t, A, B, C, kf, kb, _exp=np.exp):
    """
    Analytic solution to the dimeriztion reaction:
        A + B <=> C; (K = kf/kb)
    """
    q = -A*kf - B*kf - kb
    U = (A**2*kf**2 - 2*A*B*kf**2 + 2*A*kb*kf + B**2*kf**2 +
         2*B*kb*kf + 4*C*kb*kf + kb**2)**0.5
    expUt = _exp(U*t)
    return -(U - q)*(U + q)*(1 - 1/expUt)/(2*kf*(U + U/expUt - q + q/expUt))


def _get_Cref(t, y0, k, use_mpmath=True):
    """ convenience function for generating reference trajectory """
    if use_mpmath:
        import mpmath as mp
        mp.mp.dps = 30  # number of significant figures
        y0 = [mp.mpf(_) for _ in y0]
        k = [mp.mpf(_) for _ in k]
        _exp = np.vectorize(mp.exp)
    else:
        def _exp(x):
            return np.exp(_algebraic_sigmoid(np.asarray(x), 8, 350))
    A, B, C = y0
    kf, kb = k
    x = analytic_x(t, A, B, C, kf, kb, _exp).reshape((t.size, 1))
    dy = np.hstack((-x, -x, x))
    res = y0 + dy
    if use_mpmath:
        res = np.array(res, dtype=np.float64)
    return res


def integrate_rd(
        tend=1.9, A0=4.2, B0=3.1, C0=1.4, nt=100, t0=0.0, kf=0.9, kb=0.23,
        atol='1e-7,1e-6,1e-5', rtol='1e-6', integrator='scipy', method='bdf',
        logy=False, logt=False, num_jac=False, plot=False, savefig='None',
        splitplots=False, plotlogy=False, plotsymlogy=False, plotlogt=False,
        scale_err=1.0, scaling=1.0, verbose=False):
    """
    Runs the integration and (optionally) plots:

    - Individual concentrations as function of time
    - Reaction Quotient vs. time (with equilibrium constant as reference)
    - Numerical error commited (with tolerance span plotted)
    - Excess error committed (deviation outside tolerance span)

    Concentrations (A0, B0, C0) are taken to be in "M" (molar),
    kf in "M**-1 s**-1" and kb in "s**-1", t0 and tend in "s"
    """

    rtol = float(rtol)
    atol = list(map(float, atol.split(',')))
    if len(atol) == 1:
        atol = atol[0]
    registry = SI_base_registry.copy()
    registry['amount'] = 1.0/scaling*registry['amount']
    registry['length'] = registry['length']/10  # decimetre

    kf = kf/molar/second
    kb = kb/second

    rd = ReactionDiffusion.nondimensionalisation(
        3, [[0, 1], [2]], [[2], [0, 1]], [kf, kb], logy=logy, logt=logt,
        unit_registry=registry)

    C0 = np.array([A0, B0, C0])*molar
    if plotlogt:
        eps = 1e-16
        tout = np.logspace(np.log10(t0+eps), np.log10(tend+eps), nt)*second
    else:
        tout = np.linspace(t0, tend, nt)*second

    integr = Integration.nondimensionalisation(
        rd, C0, tout, integrator=integrator, atol=atol, rtol=rtol,
        with_jacobian=not num_jac, method=method)
    Cout = integr.with_units('Cout')
    yout, info = integr.yout, integr.info
    try:
        import mpmath
        assert mpmath  # silence pyflakes
    except ImportError:
        use_mpmath = False
    else:
        use_mpmath = True
    time_unit = get_derived_unit(registry, 'time')
    conc_unit = get_derived_unit(registry, 'concentration')
    Cref = _get_Cref(
        to_unitless(tout - tout[0], time_unit),
        to_unitless(C0, conc_unit),
        [to_unitless(kf, 1/time_unit/conc_unit),
         to_unitless(kb, 1/time_unit)],
        use_mpmath
    ).reshape((nt, 1, 3))*conc_unit
    if verbose:
        print(info)

    if plot:
        npltcols = 3 if splitplots else 1
        import matplotlib.pyplot as plt
        plt.figure(figsize=(18 if splitplots else 6, 10))

        def subplot(row=0, idx=0, adapt_yscale=True, adapt_xscale=True,
                    span_all_x=False):
            offset = idx if splitplots else 0
            ax = plt.subplot(4, 1 if span_all_x else npltcols,
                             1 + row*npltcols + offset)
            if adapt_yscale:
                if plotlogy:
                    ax.set_yscale('log')
                elif plotsymlogy:
                    ax.set_yscale('symlog')
            if adapt_xscale and plotlogt:
                ax.set_xscale('log')
            return ax

        tout_unitless = to_unitless(tout, second)
        c = 'rgb'
        for i, l in enumerate('ABC'):
            # Plot solution trajectory for i:th species
            ax_sol = subplot(0, i)
            ax_sol.plot(tout_unitless, to_unitless(Cout[:, 0, i], molar),
                        label=l, color=c[i])

            if splitplots:
                # Plot relative error
                ax_relerr = subplot(1, 1)
                ax_relerr.plot(
                    tout_unitless, Cout[:, 0, i]/Cref[:, 0, i] - 1.0,
                    label=l, color=c[i])
                ax_relerr.set_title("Relative error")
                ax_relerr.legend(loc='best', prop={'size': 11})

                # Plot absolute error
                ax_abserr = subplot(1, 2)
                ax_abserr.plot(tout_unitless, Cout[:, 0, i]-Cref[:, 0, i],
                               label=l, color=c[i])
                ax_abserr.set_title("Absolute error")
                ax_abserr.legend(loc='best', prop={'size': 11})

            # Plot absolute error
            linE = Cout[:, 0, i] - Cref[:, 0, i]
            try:
                atol_i = atol[i]
            except:
                atol_i = atol
            wtol_i = (atol_i + rtol*yout[:, 0, i])*get_derived_unit(
                rd.unit_registry, 'concentration')

            if np.any(np.abs(linE/wtol_i) > 1000):
                # Plot true curve in first plot when deviation is large enough
                # to be seen visually
                ax_sol.plot(tout_unitless, to_unitless(Cref[:, 0, i], molar),
                            label='true '+l, color=c[i], ls='--')

            ax_err = subplot(2, i)
            plot_solver_linear_error(integr, Cref, ax_err, si=i,
                                     scale_err=1/wtol_i, color=c[i], label=l)
            ax_excess = subplot(3, i, adapt_yscale=False)
            plot_solver_linear_excess_error(integr, Cref, ax_excess,
                                            si=i, color=c[i], label=l)

        # Plot Reaction Quotient vs time
        ax_q = subplot(1, span_all_x=False, adapt_yscale=False,
                       adapt_xscale=False)
        Qnum = Cout[:, 0, 2]/(Cout[:, 0, 0]*Cout[:, 0, 1])
        Qref = Cref[:, 0, 2]/(Cref[:, 0, 0]*Cref[:, 0, 1])
        ax_q.plot(tout_unitless, to_unitless(Qnum, molar**-1),
                  label='Q', color=c[i])
        if np.any(np.abs(Qnum/Qref-1) > 0.01):
            # If more than 1% error in Q, plot the reference curve too
            ax_q.plot(tout_unitless, to_unitless(Qref, molar**-1),
                      '--', label='Qref', color=c[i])
        # Plot the
        ax_q.plot((tout_unitless[0], tout_unitless[-1]),
                  [to_unitless(kf/kb, molar**-1)]*2,
                  '--k', label='K')
        ax_q.set_xlabel('t')
        ax_q.set_ylabel('[C]/([A][B]) / M**-1')
        ax_q.set_title("Transient towards equilibrium")
        ax_q.legend(loc='best', prop={'size': 11})

        for i in range(npltcols):
            subplot(0, i, adapt_yscale=False)
            plt.title('Concentration vs. time')
            plt.legend(loc='best', prop={'size': 11})
            plt.xlabel('t')
            plt.ylabel('[X]')

            subplot(2, i, adapt_yscale=False)
            plt.title('Absolute error in [{}](t) / wtol'.format('ABC'[i]))
            plt.legend(loc='best')
            plt.xlabel('t')
            ttl = '|E_i[{0}]|/(atol_i + rtol*(y0_i+yf_i)/2'
            plt.ylabel(ttl.format('ABC'[i]))
            plt.tight_layout()

            subplot(3, i, adapt_yscale=False)
            ttl = 'Excess error in [{}](t) / integrator linear error span'
            plt.title(ttl.format(
                'ABC'[i]))
            plt.legend(loc='best')
            plt.xlabel('t')
            plt.ylabel('|E_excess[{0}]| / e_span'.format('ABC'[i]))

        plt.tight_layout()
        save_and_or_show_plot(savefig=savefig)

    return yout, to_unitless(Cref, conc_unit), rd, info


if __name__ == '__main__':
    argh.dispatch_command(integrate_rd, output_file=None)
