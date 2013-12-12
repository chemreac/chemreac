#!/usr/bin/env python
# -*- coding: utf-8 -*-

#stdlib imports
import time
import json
import cPickle as pickle

#external imports
import argh
import numpy as np

#project internal imports
from reactiondiffusion.util import wrap
from reactiondiffusion.serialization import load

"""
Demo of a large chemical reaction diffusion system.
"""


def main(nt=0, tend=10.0, implementation='py', N=3,
         plot=False):
    if implementation == 'py':
        from reactiondiffusion.chem import ReactionDiffusion
    elif implementation == 'cy':
        from reactiondiffusion.cy_chem import ReactionDiffusion
    elif implementation == 'cpp':
        from reactiondiffusion.cpp_chem_wrapper import \
            PyReactionDiffusion as ReactionDiffusion
    else:
        raise NotImplementedError

    sys = load('aqueous_radiolysis'+'.json', ReactionDiffusion, N=N)


    y0_by_name = json.load(open('aqueous_radiolysis'+'.y0.json', 'rt'))
    names = json.load(open('aqueous_radiolysis'+'.names.json', 'rt'))

    # y0 with a H2 gradient
    y0 = np.array([[y0_by_name.get(k, 1e-9) if k != 'H2' else \
                    1e-3/(i+2) for k in names] for i in range(sys.N)])

    t0 = 0.0
    h = 1e-9

    if True:
        from scipy.integrate import ode
        if nt == 0: nt = 1024 # Good plotting density
        fout = np.empty(sys.n*sys.N)
        jout = np.zeros((sys.n*3+1, sys.n*sys.N), order="F")
        def f(t, y, *fargs):
            # Python function closure circumvents reallocation
            f.neval += 1
            sys.f(t, y, fout)
            return fout
        f.neval = 0
        def jac(t, y, *jargs):
            jac.neval += 1
            sys.banded_packed_jac_cmaj(t, y, jout, 1.0, False)
            return jout
        jac.neval = 0
        ode = ode(f, jac=jac)
        ode.set_integrator(
            'vode', method='bdf', atol=1e-12, rtol=1e-6,
            with_jacobian=True, lband=sys.n, uband=sys.n,
            first_step=1e-9)
        ode.set_initial_value(y0.flatten(), t0)
        tout = np.logspace(-9,np.log10(tend),nt)
        yout = np.empty((nt, sys.n*sys.N))
        texec = time.time()
        for i in range(nt):
            ode.integrate(tout[i])
            yout[i, :] = ode.y
        texec = time.time() - texec
        print("texec={0}, f.neval={1}, jac.neval={2}".format(
            texec, f.neval, jac.neval))
        xlabel = "log10(t / s)"
        ylabel = "log10(C / M)"

    else:
        from reactiondiffusion.methods import RefScipy
        m = RefScipy(sys, y0.flatten(), t0)
        texec = m.integrate(tend, nt, h=h)
        m.print_info(texec)
        tout = m.tout
        yout = m.yout[0,:,:]
        xlabel = "t / s"
        ylabel = "C / M"

    if plot:
        import matplotlib.pyplot as plt
        for subs in ('H2', 'H2O2'):
            plt.loglog(tout, yout[:,names.index(subs)],
                       label=subs+'(1)')
            plt.loglog(tout, yout[:,sys.n*(sys.N-1)+names.index(subs)],
                       label=subs+'({0})'.format(sys.N))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title("1 min post radiolysis,"+\
                  " H2 gradient ({0} bins)".format(sys.N))
        plt.legend(loc='lower left')
        plt.savefig('plot.png')#show()


argh.dispatch_command(main)
