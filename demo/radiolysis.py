#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import json
from math import log, e, exp

# external imports
import argh
import numpy as np

# project internal imports
from chemreac.serialization import load
from chemreac import ReactionDiffusion
from chemreac.integrate import integrate_scipy, integrate_sundials
from chemreac.util.plotting import plot_C_vs_t_in_bin

"""
Demo of a large chemical reaction diffusion system.
"""


def main(t0=1e-7, tend=.1, doserate=15, N=1000, nt=1024, plot=False,
         logy=False, logt=False, show=False, integrator='scipy',
         name='aqueous_radiolysis', num_jacobian=False):

    null_conc = 1e-24

    mu = 1.0  # linear attenuation
    rho = 1.0  # kg/dm3
    rd = load(name+'.json', ReactionDiffusion, N=N, logy=logy,
              logt=logt, bin_k_factor=[
                   [doserate*rho*exp(-mu*i/N)] for i in range(N)])
    y0_by_name = json.load(open(name+'.y0.json', 'rt'))

    # y0 with a H2 gradient
    y0 = np.array([[y0_by_name.get(k, null_conc) if k != 'H2' else
                    1e-3/(i+2) for k in rd.substance_names]
                   for i in range(rd.N)])

    t0 = 1e-7

    tout = np.logspace(log(t0), log(tend), nt+1, base=e)
    y = np.log(y0.flatten()) if logy else y0.flatten()
    t = np.log(tout) if logt else tout
    yout, info = {
        'scipy': integrate_scipy,
        'sundials': integrate_sundials}[integrator](
            rd, y, t, with_jacobian=(not num_jacobian))
    yout = np.exp(yout) if logy else yout

    print("texec={0}, f.neval={1}, jac.neval={2}, integrator={3}".format(
        info['texec'], info['neval_f'], info['neval_j'], integrator))

    if plot:
        import matplotlib.pyplot as plt
        bt_fmtstr = ("C(t) in bin {{0:.2g}}-{{1:.2g}} "
                     "with local doserate {}")
        ax = plt.subplot(2, 1, 1)
        plot_C_vs_t_in_bin(
            rd, tout, yout, 0, ax, substances=('H2', 'H2O2'),
            ttlfmt=bt_fmtstr.format(rd.bin_k_factor[0][0]))
        ax = plt.subplot(2, 1, 2)
        plot_C_vs_t_in_bin(
            rd, tout, yout, N-1, ax, substances=('H2', 'H2O2'),
            ttlfmt=bt_fmtstr.format(rd.bin_k_factor[N-1][0]))
        plt.tight_layout()
        if show:
            plt.show()
        else:
            plt.savefig('radiolysis.png')


if __name__ == '__main__':
    argh.dispatch_command(main)
