#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This examples plots the error in the solution
as function of number of bins. We expect different
behaviour depending on the number of stencil points used.
(N**-2, N**-4 and N**-6 for 3, 5 and 7 stencil points respectively)

"""

from __future__ import (
    print_function, division, absolute_import, unicode_literals
)

import argh
import numpy as np

from chemreac import FLAT, CYLINDRICAL, SPHERICAL, Geom_names

from analytic_diffusion import (
    flat_analytic, spherical_analytic, cylindrical_analytic, integrate_rd
)


def main(output='analytic_N_scaling.png'):
    import matplotlib.pyplot as plt
    nstencils = [3, 5, 7]
    c = 'rbk'
    m = 'osd'
    ls = ['--', ':', '-.']

    nNs = 7
    Ns = [16*2**i for i in range(nNs)]
    rates = [0, 0.1]
    for gi, geom in enumerate([FLAT, CYLINDRICAL, SPHERICAL]):
        for ri, rate in enumerate(rates):
            for si, nstencil in enumerate(nstencils):
                print(Geom_names[geom], nstencil, rate)
                tout, yout, info, rmsd_over_atol, sys = zip(*[
                    integrate_rd(N=N, nstencil=nstencil, k=rate,
                                 geom='fcs'[geom], atol=1e-8, rtol=1e-10)
                    for N in Ns])
                print('\n'.join(str(N)+': '+str(nfo) for
                                N, nfo in zip(Ns, info)))
                err = np.average(rmsd_over_atol, axis=1)
                logNs = np.log(Ns)
                logerr = np.log(err)
                p = np.polyfit(logNs[:nNs-si*2], logerr[:nNs-si*2], 1)

                plt.subplot(3, 2, gi*2 + ri + 1)
                plt.loglog(Ns, err, marker=m[si], ls='None', c=c[si])
                plt.loglog(
                    Ns[:nNs-si*2], np.exp(np.polyval(p, logNs[:nNs-si*2])),
                    ls='--', c=c[si],
                    label=str(nstencil)+': '+str(round(-p[0], 1)))
                plt.xlabel('N')
                ax = plt.gca()
                # ax.set_xticklabels(map(str, Ns))
                plt.ylabel('RMSD/atol')
                plt.legend(prop={'size': 11})
                if rate == 0:
                    plt.title('Diffusion, geom='+Geom_names[geom])
                else:
                    plt.title('Diffusion + 1 decay reaction, geom=' +
                              Geom_names[geom])

    plt.tight_layout()
    if output:
        plt.savefig(output)
    else:
        plt.show()


if __name__ == '__main__':
    argh.dispatch_command(main)
