#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Demo of non-linear fit to rate of binary reaction.
(e.g. stopped flow where one reactant is in excess)
"""

from __future__ import (absolute_import, division, print_function)

from math import ceil

import argh
import numpy as np
from numpy import exp
from scipy.optimize import curve_fit


from chemreac import ReactionDiffusion
from chemreac.integrate import run


def binary_eq_analytic(t, kf, kb, d, eps_l, c0):
    # TODO: add for c0[2] > 0
    # TODO: add for 4*A*B == (A+B+kb/kf)**2

    A, B, C = c0
    assert C == 0.0
    p = 4*A*B*kf**2
    q = kb+kf*(A+B)
    r = (q**2 - p)**0.5
    a = q+r
    b = q-r
    return eps_l*a*b*(exp(r*(t+d))-1)/(2*kf*(a*exp(r*(t+d))-b))


def fit_binary_eq(tdata, ydata, c0, **kwargs):
    def f(t, kf, kb, d, eps_l):
        return binary_eq_analytic(t, kf, kb, d, eps_l, c0)
    popt, pcov = curve_fit(f, tdata, ydata, **kwargs)
    return popt


def binary_fw_analytic(t, k, d, eps_l, c0):
    """
    k - rate coefficient
    d - delay
    eps_l - epsilon*l in Lambert-Beer's law (A=eps*b*C)
    c0 - initial conc
    """
    # dC/dt = dx/dt = k(c[0]-x)(c[1]-x) ...
    return (eps_l*(exp((c0[0]-c0[1])*k*(t+d))-1) /
            (c0[0]/c0[1]*exp((c0[0]-c0[1])*k*(t+d))-1))


def fit_binary_fw(tdata, ydata, c0, **kwargs):
    """
    Returns k, d, eps_l
    """
    def f(t, k, d, eps_l):
        return binary_fw_analytic(t, k, d, eps_l, c0)
    popt, pcov = curve_fit(f, tdata, ydata, **kwargs)
    return popt


def fit_binary_eq_rd(tdata, ydata, c0, Keq, **kwargs):
    # A + B <-> C
    rd = ReactionDiffusion(3, [[0, 1], [2]], [[2], [0, 1]], k=[0, 0])
    pconv = []

    def fit_func(tout, k_fw, tdelay, eps_l):
        pconv.append((k_fw, tdelay, eps_l))
        rd.k = [k_fw, k_fw/Keq]
        if tdelay > 0.0:
            integr = run(rd, c0, [0, tdelay])
            c1 = integr.Cout[1, 0, :]
        else:
            c1 = c0
        integr = run(rd, c1, tdelay+tout)
        return integr.Cout[:, 0, 2]*eps_l
    popt, pcov = curve_fit(fit_func, tdata, ydata, **kwargs)
    return popt, np.asarray(pconv)


def fit_binary_eq_from_temporal_abs_data(
        tdata, ydata, c0, Keq, plot_info=False,
        transient_yfrac=0.3, tail_xfrac=0.25, peak_yfrac=0.8,
        pseudo_fo=None):
    """
    Assume we follow absorbance of coloured solute B over time:

    A + B -> C        k_fw=?
    C     -> A + B    k_bw=k_fw/Keq

    Assumes C_A(t=0) > C_B(t=0) (pseudo first order guess).
    """
    concA, concB, concC = c0
    if pseudo_fo is None:
        # Could be supported by swaping concentrations...
        assert concA >= concB
        pseudo_fo = concA > concB*2

    # Guess plateau (tail-fraction)
    tf = ceil(len(ydata)*tail_xfrac)
    tf_avg = np.sum(ydata[-tf:])/tf
    eps_l_guess = tf_avg/c0[1]

    # Catch transient
    itransient = np.argwhere(ydata > ydata[0] +
                             (tf_avg-ydata[0])*peak_yfrac)[0]

    # Guess k
    Bdata = (tf_avg-ydata[:itransient])/eps_l_guess
    if pseudo_fo:
        # concA > concB
        linB = np.log(np.abs(Bdata))
        linB_lbl = "log(|[B]|)"
    else:
        # concA ~= concB
        linB = 1.0/Bdata
        linB_lbl = "1/[B]"

    lin_p = np.polyfit(tdata[:itransient], linB, 1)

    if pseudo_fo:
        k_fw_guess = -lin_p[0]
        lin_y_for_d = np.log(c0[1])
        y_for_d_lbl = 'log(B0)'
    else:
        k_fw_guess = lin_p[0]
        lin_y_for_d = 1/c0[1]
        y_for_d_lbl = '1/B0'

    d_guess = -(lin_y_for_d-lin_p[1])/lin_p[0]

    print("Psuedo first order linear fit only fw: k={}, d={}, eps_l={}".format(
        k_fw_guess, d_guess, eps_l_guess))

    k_fwnl, d_fwnl, eps_l_fwnl = fit_binary_fw(
        tdata, ydata, c0, p0=(k_fw_guess, d_guess, eps_l_guess))
    yfwnlfit = binary_fw_analytic(tdata-d_fwnl, k_fwnl, d_fwnl,
                                  eps_l_fwnl, c0)

    print("Nonlinear opt only fw: k={}, d={}, eps_l={}".format(
        k_fwnl, d_fwnl, eps_l_fwnl))

    kfw_eqnl, kbw_eqnl, d_eqnl, eps_l_eqnl = fit_binary_eq(
        tdata, ydata, c0, p0=(k_fw_guess, k_fw_guess/Keq,
                              d_guess, eps_l_guess))
    yeqnlfit = binary_eq_analytic(tdata-d_eqnl, kfw_eqnl, kbw_eqnl, d_eqnl,
                                  eps_l_eqnl, c0)

    print("Nonlinear opt equilibrium: kfw={}, kbw={}, d={}, eps_l={}".format(
        kfw_eqnl, kbw_eqnl, d_eqnl, eps_l_eqnl))

    popt, pconv = fit_binary_eq_rd(tdata, ydata, c0, Keq,
                                   p0=[kfw_eqnl, d_eqnl, eps_l_eqnl])
    print("Shooting opt fw+bw:    k={}, d={}, eps_l={}".format(
        *popt))

    if plot_info:
        import matplotlib.pyplot as plt
        # Guess of rate coefficient (k)
        plt.subplot(4, 1, 1)
        linB_t_plot = np.array([0, d_guess+tdata[itransient-1]])
        poly_linB = np.polyval((lin_p[0], lin_p[1]-d_guess*lin_p[0]),
                               linB_t_plot)
        plt.plot([tdata[0], d_guess+tdata[itransient]],
                 [lin_y_for_d, lin_y_for_d], '--', label=y_for_d_lbl)
        plt.plot(d_guess+tdata[:itransient], linB, label=linB_lbl)
        plt.plot(linB_t_plot, poly_linB,
                 label="k0={0:7.3g} d0={1:7.3g}".format(
                     k_fw_guess, d_guess))
        plt.title("Fit used for guess of rate coefficient")
        plt.legend(loc='best')

        plt.subplot(4, 1, 2)
        plt.plot(tdata, ydata, label='Input data')
        plt.plot(tdata-d_fwnl, yfwnlfit, label='Non-lin. fw fit')
        plt.legend(loc='best')

        plt.subplot(4, 1, 3)
        plt.plot(tdata, ydata, label='Input data')
        plt.plot(tdata-d_eqnl, yeqnlfit, label='Non-lin. eq fit')
        plt.legend(loc='best')

        plt.subplot(4, 1, 4)
        plt.plot(pconv[0, 0]/pconv[:, 0], label='k')
        plt.plot(pconv[0, 1]/pconv[:, 1], label='t_delay')
        plt.plot(pconv[0, 2]/pconv[:, 2], label='eps*l')
        plt.title("Convergence")
        plt.legend(loc='best')
        plt.show()

    return popt


def simulate_stopped_flow(rd, t, c0, k, noiselvl, eps_l, tdelay=None):
    if tdelay is None:
        tdelay = np.abs(np.random.normal(
            t[-1]/20, scale=t[-1]/20))
    integr = run(rd, c0, t)
    ytrue = eps_l*integr.Cout[:, 0, 2]
    skip_nt = np.argwhere(t >= tdelay)[0]
    tinp = t[:-skip_nt] if skip_nt > 0 else t
    yinp = ytrue[skip_nt:] + noiselvl*np.random.normal(
        size=len(t)-skip_nt)
    return tinp, yinp


def main(tdelay=1.0, B0=0.6, noiselvl=3e-4, nt=200, eps_l=4200.0, plot=False):
    """
    Solution:
      1. non-linear fit to:
       A) if approx equal conc A+A -> C
       B) else: A+B -> C
      2. Use as guess for guess and shoot. A + B <-> C
    """
    Keq = 10.0
    k_fw_true = 1.3
    ktrue = [1.3, k_fw_true/Keq]

    c0 = [1.0, B0, 0.0]
    ttrue = np.linspace(0, 10, nt)
    rd_eq = ReactionDiffusion(3, [[0, 1], [2]], [[2], [0, 1]], k=ktrue)
    tinp, yinp = simulate_stopped_flow(
        rd_eq, ttrue, c0, ktrue, noiselvl, eps_l, tdelay)
    k_fw_opt, d_opt, eps_l_opt = fit_binary_eq_from_temporal_abs_data(
        tinp, yinp, c0, Keq, True)

    rd_eq.k = [k_fw_opt, k_fw_opt/Keq]
    integr = run(rd_eq, c0, ttrue)
    yopt = integr.Cout[:, 0, 2]*eps_l_opt

    if plot:
        import matplotlib.pyplot as plt
        # Plot
        plt.subplot(2, 1, 1)
        plt.plot(tinp, yinp, label='Input data')
        plt.plot(ttrue-tdelay, yopt,
                 label='Shooting Opt (k={})'.format(k_fw_opt))
        plt.legend(loc='best')

        plt.subplot(2, 1, 2)
        plt.plot(tinp, yinp, label='Input data')
        # TODO: this needs to be improved...
        yquad = (B0-1/(1/B0+k_fw_opt*ttrue))*eps_l_opt
        plt.plot(ttrue-tdelay, yquad, label='Equal initial conc treatment')
        plt.legend(loc='best')
        plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
