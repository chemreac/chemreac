#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

def fit_rd(rd, tdata, ydata, c0, k_guess):
    def fit_func(tout, k):
        print k
        rd.k = [k]
        yout, info = run(rd, c0, tout)
        return yout[:,2]
    return curve_fit(fit_func, tdata, ydata, p0=k_guess)


def fit_binary_reaction_from_product(tdata, ydata, c0):
    """
    A + B -> C               k=?

    Assume C_A(t=0) > C_B(t=0) (pseduo first order guess).
    """
    from chemreac.chemistry import mk_sn_dict_from_names
    sbstncs = mk_sn_dict_from_names('ABC')
    r1 = Reaction({'A': 1, 'B': 1}, {'C': 1}, k=0.0)
    rsys = ReactionSystem([r1])
    rd = rsys.to_ReactionDiffusion(sbstncs)

    ndata = np.argwhere(ydata > ydata[-1]*0.9)[0]
    lnB = np.log(np.abs(ydata[-1]-ydata[:-ndata]))
    k_guess = -np.polyfit(tdata[:-ndata], lnB, 1)[1]
    popt, pcov = fit_rd(rd, tdata, ydata, c0, [k_guess])
    return popt


def main():
    ktrue = 0.1337
    tdata = np.linspace(0,10,1000)
    ytrue = 1-np.exp(-ktrue*tdata)
    ydata = 0.2*np.random.normal(size=len(tdata))
    c0 = [1.0, 0.1, 0.0]
    kopt = fit_binary_reaction_from_product(tdata, ydata, c0)
    print(ktrue)
    print(kopt)


if __name__ == '__main__':
    argh.dispatch_command(main)
