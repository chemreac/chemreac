#!/usr/bin/env python
# -*- coding: utf-8 -*-

#stdlib imports
import json

#external imports
import argh
import numpy as np

#project internal imports
from chemreac.serialization import load
from chemreac import ReactionDiffusion

from chemreac.integrate import run

"""
Demo of a large chemical reaction diffusion system.
"""

name = 'aqueous_radiolysis'

def main(tend=10.0, N=3, nt=1024,
         plot=False, logy=False, logt=False):

    sys = load(name+'.json', ReactionDiffusion, N=N)
    y0_by_name = json.load(open(name+'.y0.json', 'rt'))
    names = json.load(open(name+'.names.json', 'rt'))

    # y0 with a H2 gradient
    y0 = np.array([[y0_by_name.get(k, 1e-9) if k != 'H2' else \
                    1e-3/(i+2) for k in names] for i in range(sys.N)])

    t0 = 1e-7
    h = 1e-9

    tout = np.linspace(t0, tend, nt+1)

    y = np.log(y0.flatten()) if logy else y0.flatten()
    t = np.log(tout) if logt else tout
    yout, info = run(sys, y, t)
    if logy: yout = np.exp(yout)

    print("texec={0}, f.neval={1}, jac.neval={2}".format(
        info['texec'], info['neval_f'], info['neval_j']))

    if plot:
        import matplotlib.pyplot as plt
        if logy:
            if logt:
                plt_cb = plt.plot
            else:
                plt_cb = plt.semilogy
        else:
            if logt:
                plt_cb = plt.semilogx
            else:
                plt_cb = plt.loglog

        for subs in ('H2', 'H2O2'):
            plt_cb(t, yout[:,names.index(subs)],
                       label=subs+'(1)')
            plt_cb(t, yout[:,sys.n*(sys.N-1)+names.index(subs)],
                       label=subs+'({0})'.format(sys.N))
        xlabel = "log10(t / s)"
        ylabel = "log10(C / M)"
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title("1 min post radiolysis,"+\
                  " H2 gradient ({0} bins)".format(sys.N))
        plt.legend(loc='lower left')
        plt.savefig('plot.png')#show()


argh.dispatch_command(main)
