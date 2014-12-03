#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Four species two reactions
--------------------------

:download:`examples/four_species.py` demonstrates how
to use the plotting utilities ``plot_per_reaction_contribution``
and ``plot_jacobian``. The reaction system is as follows:

.. math::

    A &\rightarrow B ~~~ & k_1=0.05 \\
    2C + B &\rightarrow D + B ~~~ & k_2=3.0


::

 $ python four_species.py --help

.. exec::
   echo "::\\n\\n"
   python examples/examples/four_species.py --help | sed "s/^/   /"


Here is an example generated by:

::

 $ python four_species.py --plot --savefig four_species.png

.. image:: ../_generated/four_species.png

Jacobian:

.. image:: ../_generated/four_species_jacobian.png

Per reaction contribution:

.. image:: ../_generated/four_species_per_reaction.png

"""

import os

import argh
import numpy as np

from chemreac import DENSE, BANDED, SPARSE
from chemreac.integrate import run
from chemreac.serialization import load
from chemreac.util.plotting import (
    coloured_spy, plot_jacobian, plot_per_reaction_contribution,
    save_and_or_show_plot
)

# A      -> B          k1=0.05
# 2C + B -> D + B      k2=3.0


def integrate_rd(tend=10.0, N=1, nt=500, jac_spy=False, mode=None,
                 logy=False, logt=False, plot=False, savefig='None',
                 verbose=False):
    """
    Integrates the reaction system defined by
    :download:`four_species.json <examples/four_species.json>`
    """
    rd = load(os.path.join(os.path.dirname(
        __file__), 'four_species.json'), N=N, x=N, logy=logy, logt=logt)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    y0 = np.concatenate([y0/(i+1)*(0.25*i**2+1) for i in range(N)])
    t0 = 1e-10

    if mode is None:
        if rd.N == 1:
            mode = DENSE
        elif rd.N > 1:
            mode = BANDED
    else:
        mode = int(mode)

    import matplotlib.pyplot as plt
    if jac_spy:
        fout = np.empty(rd.n*rd.N)
        rd.f(t0, y0, fout)
        print(fout)
        if mode == DENSE:
            jout = np.zeros((rd.n*rd.N, rd.n*rd.N), order='F')
            rd.dense_jac_cmaj(t0, y0, jout)
            coloured_spy(np.log(np.abs(jout)))
        elif mode == BANDED:
            # note rd.n*3 needed in call from scipy.integrate.ode
            jout = np.zeros((rd.n*2+1, rd.n*rd.N), order='F')
            rd.banded_packed_jac_cmaj(t0, y0, jout)
            coloured_spy(np.log(np.abs(jout)))
        print(jout)
        plt.show()
    else:
        tout = np.linspace(t0, tend, nt)
        integr = run(rd, y0, tout)
        Cout = integr.Cout
        if verbose:
            print(integr.info)
        if plot:
            plt.figure(figsize=(6, 4))
            for i, l in enumerate('ABCD'):
                plt.plot(tout, Cout[:, 0, i], label=l)
            plt.title("Time evolution of concentrations")
            plt.legend()
            save_and_or_show_plot(savefig=savefig)

            plt.figure(figsize=(6, 10))
            plot_jacobian(
                rd,
                np.log(tout) if rd.logt else tout,
                np.log(Cout) if rd.logy else Cout,
                'ABCD',
                lintreshy=1e-10
            )
            plt.tight_layout()
            if savefig != 'None':
                base, ext = os.path.splitext(savefig)
                savefig = base + '_jacobian' + ext
            save_and_or_show_plot(savefig=savefig)

            plt.figure(figsize=(6, 10))
            plot_per_reaction_contribution(
                rd,
                np.log(tout) if rd.logt else tout,
                np.log(Cout) if rd.logy else Cout,
                'ABCD'
            )
            plt.tight_layout()
            if savefig != 'None':
                savefig = base + '_per_reaction' + ext
            save_and_or_show_plot(savefig=savefig)


if __name__ == '__main__':
    argh.dispatch_command(integrate_rd, output_file=None)
