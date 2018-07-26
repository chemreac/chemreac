#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

import argh
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from decay import get_Cref


def main(tend=2.0, A0=1.0, nt=67, t0=0.0,
         rates='3.40715,4.0'):
    k = list(map(float, rates.split(',')))
    n = len(k)+1
    if n > 4:
        raise ValueError("Max 3 consequtive decays supported at the moment.")
    tout = np.linspace(t0, tend, nt)
    y0 = np.zeros(n)
    y0[0] = A0
    Cref = get_Cref(k, y0, tout - tout[0]).reshape((nt, 1, n))

    plt.xkcd()
    fig = plt.figure(figsize=(2, 2), dpi=100)
    ax = plt.subplot(1, 1, 1)
    for i, l in enumerate('ABC'[:n]):
        ax.plot(tout, Cref[:, 0, i], label=l, color='rbg'[i])
    ax.xaxis.set_tick_params(width=1)
    ax.yaxis.set_tick_params(width=1)
    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 0.5, 1])
    ax.set_facecolor((0.9, 0.9, 0.9))
    fig.patch.set_facecolor((1.0, 1.0, 1.0, 0.0))
    #ax.text(.35, 0.5, r'A $\rightarrow$ B $\rightarrow$ C', fontsize=9)
    plt.title('chemreac', fontsize=21)
    plt.legend(loc='best', prop={'size': 10})
    plt.tight_layout()
    plt.savefig('chemreac_logo.svg', transparent=True)
    plt.savefig('chemreac_logo.png', transparent=True)


if __name__ == '__main__':
    argh.dispatch_command(main)
