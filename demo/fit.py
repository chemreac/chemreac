#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac.integrate import run
from chemreac.chemistry import Reaction, ReactionSystem

"""
Demo of non-linear fit to rate of binary reaction.
(e.g. stopped flow where one reactant is in excess)
"""


import numpy as np
from scipy.optimize import curve_fit


def fit_rd(rd, tdata, ydata, c0, param_guess):
    pconv = []
    def fit_func(tout, k, tdelay):
        pconv.append((k, tdelay))
        rd.k = [k]
        if tdelay > 0.0:
            c1, info = run(rd, c0, [0, tdelay])
            c1 = c1[1,:]
        else:
            c1 = c0
        yout, info = run(rd, c1, tdelay+tout)
        return yout[:,2]
    popt, pcov = curve_fit(fit_func, tdata, ydata, p0=param_guess)
    return popt, np.asarray(pconv)


def _get_rd():
    from chemreac.chemistry import mk_sn_dict_from_names
    sbstncs = mk_sn_dict_from_names('ABC')
    r1 = Reaction({'A': 1, 'B': 1}, {'C': 1}, k=0.0)
    rsys = ReactionSystem([r1])
    return rsys.to_ReactionDiffusion(sbstncs)


def fit_binary_reaction_from_product(
        tdata, ydata, c0, plot_info=False, 
        transient_yfrac=0.3, tail_xrfrac=4, peak_yfrac=0.8):
    """
    A + B -> C        k=?

    Assumes C_A(t=0) > C_B(t=0) (pseduo first order guess).
    """
    concA, concB, concC = c0
    assert concA >= concB # Could be supported by swaping concentrations...
    pseudo_fo = concA > concB*2
    
    # Guess plateau
    tf = len(ydata)//tail_xrfrac
    tf_avg = np.sum(ydata[-tf:])/tf

    # Guess delay
    hf = np.argwhere(ydata > ydata[0]+(tf_avg-ydata[0])*transient_yfrac)[0]
    d_p = np.polyfit(tdata[:hf], ydata[:hf], 1)
    d_guess = ydata[0]/d_p[0]

    # Catch transient
    itransient = np.argwhere(ydata > tf_avg*peak_yfrac)[0]

    # Guess k
    Bdata = tf_avg-ydata[:itransient]
    if pseudo_fo:
        # concA > concB
        linB = np.log(np.abs(Bdata))
        linB_lbl = "log(|[B]|)"
    else:
        # concA ~= concB
        linB = 1.0/Bdata
        linB_lbl = "1/[B]"


    k_p = np.polyfit(tdata[:itransient], linB, 1)

    if pseudo_fo:
        k_guess = -k_p[0]
    else:
        k_guess = k_p[0]

    rd = _get_rd()
    popt, pconv = fit_rd(rd, tdata, ydata, c0, [k_guess, d_guess])

    if plot_info:
        # Guess of delay
        plt.subplot(3, 1, 1)
        plt.plot(d_guess+tdata[:hf], ydata[:hf], 'sr')
        plt.plot(d_guess+tdata[hf:], ydata[hf:], 'sk')
        nt_plot=min(hf*3, len(ydata)//2)
        d_t_plot = np.array([0, d_guess+tdata[nt_plot-1]])
        plt.plot(d_t_plot, np.polyval((d_p[0], d_p[1]-d_guess*d_p[0]), d_t_plot), '-r', 
                 label="delay_guess={}".format(d_guess))
        plt.title("Fit used for guess of delay")
        plt.ylim((0, plt.ylim()[1]))
        plt.legend(loc='best')

        # Guess of rate coefficient (k)
        plt.subplot(3, 1, 2)
        linB_t_plot = np.array([0, d_guess+tdata[itransient-1]])
        poly_linB = np.polyval((k_p[0], k_p[1]-d_guess*k_p[0]), linB_t_plot)
        plt.plot(d_guess+tdata[:itransient], linB, label=linB_lbl)
        plt.plot(linB_t_plot, poly_linB, label="k_guess={}".format(k_guess))
        plt.title("Fit used for guess of rate coefficient")
        plt.legend(loc='best')

        plt.subplot(3, 1, 3)
        plt.plot(pconv[:,0], label='k')
        plt.plot(pconv[:,1], label='t_delay')
        plt.title("Convergence")
        plt.legend(loc='best')
        plt.show()

    return popt


def y_for_k(tdata, c0, k):
    rd = _get_rd()
    rd.k = [k]
    yout, info = run(rd, c0, tdata)
    return yout[:,2]


def main():
    ktrue = 0.5
    nt = 200
    ttrue = np.linspace(0,10,nt)
    c0 = [1.0, 1.0-1e-9, 0.0]
    ytrue = y_for_k(ttrue, c0, ktrue)
    ypseudo = c0[1]*(1-np.exp(-ktrue*ttrue))

    # delay before meassurement
    tdelay = np.abs(np.random.normal(1.0))
    skip_nt = np.argwhere(ttrue > tdelay)[0]
    yinp = ytrue[skip_nt:] + 0.0000003*np.random.normal(size=len(ttrue)-skip_nt)
    kopt, dopt = fit_binary_reaction_from_product(
        ttrue[:-skip_nt], yinp, c0, True)
    yopt = y_for_k(ttrue, c0, kopt)

    # Plot
    plt.plot(ttrue[:-skip_nt], yinp, label='Input data')
    plt.plot(ttrue-tdelay, ypseudo, label='Pseudo-first order treatment')
    plt.plot(ttrue-tdelay, yopt, label='Opt (k={})'.format(kopt))
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
