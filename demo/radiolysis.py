#!/usr/bin/env python
# -*- coding: utf-8 -*-

#stdlib imports
import json
from math import log, e

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

def main(tend=10.0, N=2, nt=1024,
         plot=False, logy=False, logt=False, show=False):

    sys = load(name+'.json', ReactionDiffusion, N=N, logy=logy, logt=logt)
    y0_by_name = json.load(open(name+'.y0.json', 'rt'))
    names = json.load(open(name+'.names.json', 'rt'))

    # y0 with a H2 gradient
    y0 = np.array([[y0_by_name.get(k, 1e-9) if k != 'H2' else \
                    1e-3/(i+2) for k in names] for i in range(sys.N)])

    t0 = 1e-7

    tout = np.logspace(log(t0), log(tend), nt+1, base=e)
    y = np.log(y0.flatten()) if logy else y0.flatten()
    t = np.log(tout) if logt else tout
    yout, info = run(sys, y, t)
    if logy: yout = np.exp(yout)
    #if logt: tout = np.exp(tout)

    print("texec={0}, f.neval={1}, jac.neval={2}".format(
        info['texec'], info['neval_f'], info['neval_j']))

    if plot:
        import matplotlib.pyplot as plt
        for subs in ('H2', 'H2O2'):
            plt.loglog(tout, yout[:,names.index(subs)],
                       label=subs+'(1)')
            plt.loglog(tout, yout[:,sys.n*(sys.N-1)+names.index(subs)],
                       label=subs+'({0})'.format(sys.N))
        xlabel = "t / s"
        ylabel = "C / M"
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(("{} s post radiolysis,"+\
                   " H2 gradient ({} bins)").format(tend, sys.N))
        plt.legend(loc='best')
        if show:
            plt.show()
        else:
            plt.savefig('radiolysis.png')


argh.dispatch_command(main)
