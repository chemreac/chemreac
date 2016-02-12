#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The fruit fly of stiff numerical chemical kinetics problems.

$ python robertson.py -A 1.0 -B 1e-20 -C 1e-20 --t0 0 --plot --tend 3e10 --nt \
1024 --logt --logy --verbose

"""
from __future__ import (absolute_import, division, print_function)

from math import log10
import numpy as np
from chemreac import ReactionDiffusion
from chemreac.chemistry import Reaction, ReactionSystem
from chemreac.integrate import run
from chemreac.util.analysis import suggest_t0
from chemreac.util.plotting import (
    plot_C_vs_t, save_and_or_show_plot, plot_faded_time
)


def get_reactions(rates):
    """
        A -> B
    B + C -> A + C
    B + B -> C
    """
    return (
        Reaction({'A': 1}, {'B': 1}, rates[0]),
        Reaction({'B': 1, 'C': 1}, {'A': 1, 'C': 1}, rates[1]),
        Reaction({'B': 2}, {'B': 1, 'C': 1}, rates[2])
    )


def integrate_rd(tend=1e2, A0=1.0, B0=0.0, C0=0.0, k1=0.04, k2=1e4, k3=3e7,
                 t0=1e2, nt=100, N=1, nstencil=3, logt=False, logy=False,
                 plot=False, savefig='None', verbose=False, dump_expr='False',
                 use_chempy=False, D=2e-3):
    if N == 1:
        init_conc = (A0, B0, C0)
    else:
        init_conc = np.tile((A0, B0, C0), (N, 1))
        init_conc /= np.linspace(1, 2, N).reshape((N, 1))**.5

    rsys = ReactionSystem(get_reactions((k1, k2, k3)), 'ABC')
    if verbose:
        print([str(_) for _ in rsys.rxns])
    if use_chempy:
        from chempy.kinetics.ode import get_odesys
        odesys = get_odesys(rsys, include_params=True)
        if N != 1:
            raise ValueError("ChemPy does not support diffusion")
        odesys.integrate(np.logspace(log10(t0), log10(tend)), init_conc)
        if plot:
            odesys.plot_result(xscale='log', yscale='log')
        result = None
    else:
        rd = ReactionDiffusion.from_ReactionSystem(
            rsys, N=N, nstencil=1 if N == 1 else nstencil, logt=logt,
            logy=logy, D=[D/2, D/3, D/5])

        if dump_expr.lower() not in ('false', '0'):
            from chemreac.symbolic import SymRD
            import sympy as sp
            cb = {'latex': sp.latex,
                  'ccode': sp.ccode}.get(dump_expr.lower(), str)
            srd = SymRD.from_rd(rd, k=sp.symbols('k:3'))
            print('dydx:')
            print('\n'.join(map(cb, srd._f)))
            print('jac:')
            for ri, row in enumerate(srd.jacobian.tolist()):
                for ci, expr in enumerate(row):
                    if expr == 0:
                        continue
                    print(ri, ci, cb(expr))
            return None

        if t0 == 0 and logt:
            t0 = 1e-3*suggest_t0(rd, init_conc)
            if verbose:
                print("Using t0 = %12.5g" % t0)

        t = np.logspace(np.log10(t0), np.log10(tend), nt)
        print(t[0], t[-1])
        integr = run(rd, init_conc, t)

        if verbose:
            import pprint
            pprint.pprint(integr.info)

        if plot:
            if N == 1:
                plot_C_vs_t(integr, xscale='log', yscale='log')
            else:
                import matplotlib.pyplot as plt
                for idx, name in enumerate('ABC', 1):
                    plt.subplot(1, 3, idx)
                    rgb = [.5, .5, .5]
                    rgb[idx-1] = 1
                    plot_faded_time(integr, name, rgb=rgb, log_color=True)
        result = integr

    if plot:
        save_and_or_show_plot(savefig=savefig)

    return result


if __name__ == '__main__':
    import argh
    argh.dispatch_command(integrate_rd)
