#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The fruit fly of stiff numerical chemical kinetics problems.

$ python robertson.py -A 1.0 -B 1e-20 -C 1e-20 --t0 0 --plot --tend 3e10 --nt \
1024 --logt --logy --verbose
"""

import numpy as np
from chemreac.chemistry import mk_sn_dict_from_names, Reaction, ReactionSystem
from chemreac.integrate import run
from chemreac.util.analysis import suggest_t0
from chemreac.util.plotting import plot_C_vs_t, save_and_or_show_plot


def get_reaction_system(substances, rates):
    """
        A -> B
    B + C -> A + C
    B + B -> C
    """
    A, B, C = substances
    reactions = (Reaction({A: 1}, {B: 1}, k=rates[0]),
                 Reaction({B: 1, C: 1}, {A: 1, C: 1}, k=rates[1]),
                 Reaction({B: 2}, {C: 1}, k=rates[2]))
    return ReactionSystem(reactions)


def integrate_rd(tend=7200.0, A0=1.0, B0=0.0, C0=0.0, k1=0.04, k2=1e4, k3=3e7,
                 t0=0.0, nt=100, logt=False, logy=False, plot=False,
                 savefig='None', verbose=False, dump_latex=False):
    sn_dict = mk_sn_dict_from_names('ABC', mass=(1, 1, 2))
    k = (k1, k2, k3)
    init_conc = (A0, B0, C0)
    reaction_system = get_reaction_system(sn_dict.values(), k)
    print([str(_) for _ in reaction_system.rxns])
    rd = reaction_system.to_ReactionDiffusion(logt=logt, logy=logy)
    if dump_latex:
        from chemreac.symbolic import SymRD
        import sympy as sp
        srd = SymRD.from_rd(rd)
        print('dydx:')
        print('\n'.join(map(sp.latex, srd._f)))
        print('jac:')
        for ri, row in enumerate(srd.jacobian.tolist()):
            for ci, expr in enumerate(row):
                if expr == 0:
                    continue
                print(ri, ci, sp.latex(expr))
        return None
    if t0 == 0 and logt:
        t0 = 1e-3*suggest_t0(rd, init_conc)
        if verbose:
            print("Using t0 = %12.5g" % t0)
    t = np.logspace(np.log10(t0), np.log10(tend), nt)
    integr = run(rd, init_conc, t)
    if verbose:
        print(integr.info)
    if plot:
        plot_C_vs_t(integr, labels=sn_dict.keys(),
                    xscale='log', yscale='log')
        save_and_or_show_plot(savefig=savefig)
    return integr


if __name__ == '__main__':
    import argh
    argh.dispatch_command(integrate_rd)
