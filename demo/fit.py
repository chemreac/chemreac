#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac.integrate import run
from chemreac.chemistry import Reaction, ReactionSystem

"""
Demo of non-linear fit to
"""


import numpy as np
from scipy.optimize import curve_fit

def fit_rd(rd, tdata, ydata, c0, param_guess):
    def fit_func(tout, k, tdelay):
        print k
        rd.k = [k]
        yout, info = run(rd, c0, tdelay+tout)
        return yout[:,2]
    return curve_fit(fit_func, tdata, ydata, p0=param_guess)


def _get_rd():
    from chemreac.chemistry import mk_sn_dict_from_names
    sbstncs = mk_sn_dict_from_names('ABC')
    r1 = Reaction({'A': 1, 'B': 1}, {'C': 1}, k=0.0)
    rsys = ReactionSystem([r1])
    return rsys.to_ReactionDiffusion(sbstncs)


def fit_binary_reaction_from_product(
        tdata, ydata, c0, plot_pseudo_guess=False, 
        head_frac=10, tail_frac=4, peak_frac=0.8):
    """
    A + B -> C               k=?

    Assume C_A(t=0) > C_B(t=0) (pseduo first order guess).
    """

    # Guess delay
    hf = len(ydata)//head_frac
    p = np.polyfit(tdata[:hf], ydata[:hf], 1)
    d_guess = ydata[0]/p[0]

    # Guess plateau
    tf = len(ydata)//tail_frac
    tf_avg = np.sum(ydata[-tf:])/tf

    # Catch transient
    ndata = np.argwhere(ydata > tf_avg*peak_frac)[0]
    lnB = np.log(np.abs(tf_avg-ydata[:-ndata]))

    # Guess k
    p = np.polyfit(tdata[:-ndata], lnB, 1)
    k_guess = -p[0]


    if plot_pseudo_guess:
        poly_lnB = np.polyval(p, tdata[:-ndata])
        plt.plot(tdata[:-ndata], lnB)
        plt.plot(tdata[:-ndata], poly_lnB, label="k={}".format(k_guess))
        plt.legend()
        plt.show()

    rd = _get_rd()
    popt, pcov = fit_rd(rd, tdata, ydata, c0, [k_guess, d_guess])
    return popt

def y_for_k(tdata, c0, k):
    rd = _get_rd()
    rd.k = [k]
    yout, info = run(rd, c0, tdata)
    return yout[:,2]

def main():
    ktrue = 0.5
    nt = 1000
    ttrue = np.linspace(0,10,1000)
    c0 = [1.0, 0.1, 0.0]
    ytrue = y_for_k(ttrue, c0, ktrue)
    ypseudo = c0[1]*(1-np.exp(-ktrue*ttrue))

    # delay before meassurement
    tdelay = np.abs(np.random.normal(1.0))
    skip_nt = np.argwhere(ttrue > tdelay)[0]
    yinp = ytrue[skip_nt:] + 0.003*np.random.normal(size=len(ttrue)-skip_nt)
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
