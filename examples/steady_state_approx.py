#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Steady state approximation
--------------------------

:download:`examples/steady_state_approx.py` shows how you can estimate
errors commited when assuming steady state for simple systems.
We will consider the following system:

::

    A     ->  B      k1
        B ->  C      k2
    A + B ->  B + C  k3

    dA/dt = -k1*A          - k3*A*B
    dB/dt =  k1*A - k2*B
    dC/dt =         k2*B   + k3*A*B


The rate expressions are from mass action and hence
we are conserving mass:

.. math ::

    A+B+C = A_0 + B_0 + C_0


sum of derivatives = 0 (already satisfied)


For initial concentrations of A much larger than B we have:


Steady state assumption for B (A0 >> B0):

.. math ::

    B &= \frac{k_1 A}{k_2} \\
    \frac{dA}{dt} &= -k_1\left(A + \frac{k_3}{k_2}A^2 \right) \\
    \log{A} - \log{\left( \frac{k_3}{k_2}A + 1 \right)} &= -k_1t +
        \left(\log{A_0} - \log{\left( \frac{k_3}{k_2}A_0 + 1\right)}\right) \\
    \frac{A}{\frac{k_3}{k_2}A + 1} &= f(t)


using

.. math ::

    f(t) = \frac{A_0}{\frac{k_3}{k_2}A_0 + 1}e^{-k_1 t}


we get

.. math ::

    A(1 - f(t)*\frac{k_3}{k_2}) &= f(t) \\
    A &= \frac{f(t)}{1 - f(t)\frac{k_3}{k_2}} [*] \\
    A &= \frac{1}{\frac{1}{f(t)} - \frac{k_3}{k_2}} \\
    A &= \frac{1}{\frac{\frac{k_3}{k_2}A_0 + 1}{A_0}e^{k_1 t} -
        \frac{k_3}{k_2}} \\


where we note:

.. math ::

    [*] if \frac{k_3}{k_2} << 1: \\
    A = f(t)


"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

import argh
import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import run
from chemreac.util.plotting import plot_C_vs_t_in_bin, save_and_or_show_plot


def integrate_rd(tend=1.0, k1=7e-1, k2=3e2, k3=7.0,
                 A0=1.0, B0=0.0, C0=0.0, nt=512,
                 plot=False, savefig='None'):
    """
    Runs integration and (optionally) generates plots.


    Examples
    --------
    ::

       $ python steady_state_approx.py --plot --savefig steady_state_approx.png


    .. image:: ../_generated/steady_state_approx.png


    ::

       $ python steady_state_approx.py --plot --savefig \
           steady_state_approx.html


    :download:`../_generated/steady_state_approx.html`

    """
    def f(t):
        return A0/(k3/k2*A0 + 1)*np.exp(-k1*t)

    y = [A0, B0, C0]
    k = [k1, k2, k3]
    rd = ReactionDiffusion(
        3, [[0], [1], [0, 1]], [[1], [2], [1, 2]], k)
    t = np.linspace(0, tend, nt)
    yout, info = run(rd, y, t)
    A_ssB = 1/(1/f(t) - k3/k2)
    A_ssB_2fast = f(t)

    if plot:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(6, 10))

        ax = plt.subplot(3, 1, 1)
        plot_C_vs_t_in_bin(rd, t, yout, ax=ax)
        plt.subplot(3, 1, 2)
        plt.plot(t, yout[:, 0, 0] - A_ssB,
                 label="Abs. err. in A assuming steady state of B")
        plt.legend(loc='best', prop={'size': 11})
        plt.subplot(3, 1, 3)
        plt.plot(t, yout[:, 0, 0] - A_ssB_2fast,
                 label="Abs. err. in A when also assuming k2 >> k3")
        plt.legend(loc='best', prop={'size': 11})
        plt.tight_layout()

        save_and_or_show_plot(savefig=savefig)

    ydot = lambda x, y: (-k1*y[0] - k3*y[0]*y[1], k1*y[0] - k2*y[1],
                         k2*y[1] + k3*y[0]*y[1])
    return t, yout, A_ssB, A_ssB_2fast, ydot


if __name__ == '__main__':
    argh.dispatch_command(integrate_rd, output_file=None)
