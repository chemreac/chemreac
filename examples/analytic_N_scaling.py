#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Analytic error scaling vs. number of bins
-----------------------------------------

:download:`examples/analytic_N_scaling.py` plots the error in the solution
as function of number of bins. We expect different
behaviour depending on the number of stencil points used.
(N**-2, N**-4 and N**-6 for 3, 5 and 7 stencil points respectively)

::

 $ python analytic_N_scaling.py --help

.. exec::
   echo "::\\n\\n"
   python examples/examples/analytic_N_scaling.py --help | sed "s/^/   /"


Here is an example generated by:

::

 $ python analytic_N_scaling.py --nNs 6 --plot --savefig analytic_N_scaling.png


.. image:: ../_generated/analytic_N_scaling.png

"""

from __future__ import print_function, division, absolute_import

from itertools import product
from collections import defaultdict

import argh
import numpy as np

from chemreac import FLAT, CYLINDRICAL, SPHERICAL, Geom_names
from chemreac.util.plotting import save_and_or_show_plot

from analytic_diffusion import (
    integrate_rd
)


def main(plot=False, savefig='None', nNs=7, Ns=None, nspecies='1,2,3',
         nfit='7,5,4', keys='nfev,njev,texec', ylims='None'):
    import matplotlib.pyplot as plt
    nfit = [float(_) for _ in nfit.split(',')]
    nstencils = range(3, 2 + 2*len(nfit), 2)
    c = 'rbkg'
    m = 'osdx'
    keys = keys.split(',')
    if Ns is None:
        Ns = [10*(2**i) for i in range(nNs)]
    else:
        Ns = list(map(int, Ns.split(',')))
        nNs = len(Ns)

    nspecies = list(map(int, nspecies.split(',')))
    if plot:
        figs = [plt.figure(idx, figsize=(4*len(nspecies), 10))
                for idx in range(len(keys)+1)]
    geoms = [FLAT, CYLINDRICAL, SPHERICAL]
    for gi, geom in enumerate(geoms):
        for ni, ns in enumerate(nspecies):
            def _ttl_labels(ylbl):
                plt.ylabel(ylbl)
                plt.xlabel('N')
                if ns == 1:
                    plt.title('{}'.format(Geom_names[geom]))
                else:
                    plt.title('{}, {} decay(s)'.format(
                              Geom_names[geom], ns - 1))

            for si, nstencil in enumerate(nstencils):
                print(Geom_names[geom], nstencil, ns)
                tout, yout, info, rmsd_over_atol, sys = zip(*[
                    integrate_rd(N=N, nstencil=nstencil, nspecies=ns,
                                 geom='fcs'[geom], atol=1e-8, rtol=1e-10)
                    for N in Ns])
                print('\n'.join(str(N)+': '+str(nfo) for
                                N, nfo in zip(Ns, info)))
                err = np.average(rmsd_over_atol, axis=1)
                logNs = np.log(Ns)
                logerr = np.log(err)

                if plot:
                    p = np.polyfit(logNs[:nfit[si]], logerr[:nfit[si]], 1)
                    plt.figure(0)
                    ax = plt.subplot(len(geoms), len(nspecies),
                                     gi*len(nspecies) + ni + 1)
                    ax.set_xscale('log', basex=2)
                    ax.set_yscale('log', basey=2)
                    ax.plot(Ns, err, marker=m[si], ls='None', c=c[si])
                    ax.plot(
                        Ns[:nNs-si], np.exp(np.polyval(p, logNs[:nNs-si])),
                        ls='--', c=c[si],
                        label=str(nstencil)+': '+str(round(-p[0], 1)))
                    ax = plt.gca()
                    # ax.set_xticklabels(map(str, Ns))
                    plt.legend(loc='upper right', prop={'size': 10})
                    _ttl_labels('RMSD/atol')

                    for idx, key in enumerate(keys, 1):
                        plt.figure(idx)
                        ax = plt.subplot(len(geoms), len(nspecies),
                                         gi*len(nspecies) + ni + 1)
                        ax.set_xscale('log', basex=2)
                        ax.set_yscale('log')  # , basey=10)
                        ax.plot(Ns, [nfo[key] for nfo in info], marker=m[si],
                                ls='None', c=c[si], label=str(nstencil))
                        plt.legend(loc='upper left', prop={'size': 10},
                                   numpoints=1)
                        _ttl_labels(key)

    if plot:
        # Shared ylim
        gini = list(product(range(len(geoms)), range(len(nspecies))))
        if ylims == 'None':
            ylims = defaultdict(lambda: [float('nan'), float('nan')])
            for fi in range(len(keys)+1):
                plt.figure(fi)
                for gi, ni in gini:
                    ax = plt.subplot(len(geoms), len(nspecies),
                                     gi*len(nspecies) + ni + 1)
                    ylim = ax.get_ylim()
                    ylims[fi][0] = min(ylim[0], ylims[fi][0])
                    ylims[fi][1] = max(ylim[1], ylims[fi][1])
        else:
            ylims = [[float(_) for _ in ylim.split(',')]
                     for ylim in ylims.split(';')]
        for fi in range(len(keys)+1):
            plt.figure(fi)
            for gi, ni in gini:
                ax = plt.subplot(len(geoms), len(nspecies),
                                 gi*len(nspecies) + ni + 1)
                ax.set_ylim(ylims[fi])

        # Save figs
        plt.figure(0)
        plt.tight_layout()
        save_and_or_show_plot(savefig=savefig)
        for fi, key in enumerate(keys, 1):
            plt.figure(fi)
            plt.tight_layout()
            save_and_or_show_plot(savefig='.'.join([
                ('' if i > 0 else key + '_') + p for i, p
                in enumerate(savefig.split('.'))
            ]))


if __name__ == '__main__':
    argh.dispatch_command(main, output_file=None)
