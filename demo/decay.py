#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac import ReactionDiffusion
from chemreac.integrate import run


def integrate_rd(tend=2.0, A0=3.14, nt=67, small=20, logy=False, logt=False,
         plot=False, atolA=1e-6, atolB=1e-6, rtolA=1e-6, rtolB=1e-6):
    """ A -> B """
    rd = ReactionDiffusion(2, [[0]], [[1]], [1.0], 1, logy=logy, logt=logt)
    B0 = 10**(-small)
    y0 = np.array([A0, B0])
    t0 = 1e-10
    tout = np.linspace(t0, tend, nt)
    yref = np.column_stack((y0[0]*np.exp(-tout), y0[1]+y0[0]*(1-np.exp(-tout))))

    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(rd, y, t, atol=[atolA, atolB], rtol=[rtolA, rtolB])
    yout = np.exp(yout) if logy else yout
    yout = yout[:,0,:]

    if plot:
        for i, l in enumerate('AB'):
            plt.subplot(2, 1, 1)
            plt.plot(tout, yout[:, i], label=l)
            plt.subplot(2, 1, 2)
            try:
                atol = info['atol'][i]
            except:
                atol = info['atol']
            plt.plot(tout, (yout[:, i]-yref[:, i])/atol, label=l)

        plt.subplot(2, 1, 1)
        plt.title('C(t)')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.ylabel('C')
        plt.subplot(2, 1, 2)
        plt.title('Error in C(t)')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.ylabel('abs. error in C / atol')
        plt.tight_layout()
        plt.show()

    return yout, yref, rd, info

if __name__ == '__main__':
    argh.dispatch_command(integrate_rd)
