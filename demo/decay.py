#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import argh
import numpy as np
import matplotlib.pyplot as plt

from chemreac import ReactionDiffusion
from chemreac.integrate import run


def main(tend=2.0, nt=500, logy=False, logt=False):
    """ A -> B """
    sys = ReactionDiffusion(2, [[0]], [[1]], [1.0], 1, logy=logy, logt=logt)

    y0 = np.array([1.0, 1e-24])
    t0 = 1e-10
    tout = np.linspace(t0, tend, nt)
    yref = np.column_stack((np.exp(-tout), 1-np.exp(-tout)))

    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(sys, y, t)
    if logy: yout = np.exp(yout)
    for i,l in enumerate('AB'):
        plt.subplot(2,1,1)
        plt.plot(tout, yout[:,i], label=l)
        plt.subplot(2,1,2)
        plt.plot(tout, yout[:,i]-yref[:,i], label=l)

    plt.subplot(2,1,1)
    plt.title('C(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.ylabel('C')
    plt.subplot(2,1,2)
    plt.title('Error in C(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.ylabel('abs. error in C')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    argh.dispatch_command(main)
