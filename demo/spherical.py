#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argh
import numpy as np

from chemreac.cpp_chem_wrapper import \
    PyReactionDiffusion as ReactionDiffusion

from chemreac import BANDED
from chemreac.integrate import run

"""
Demo of chemical reaction diffusion system.
"""

# <spherical.png>

def main(tend=10.0, N=1, nt=50, show=False):
    sys = ReactionDiffusion(1, []N=N)

    tout, yout, info = run(sys, y0, t0, tend, nt)
    if plot:
        for i,l in enumerate('ABCD'):
            plt.plot(tout, yout[:,i], label=l)
        if show:
            plt.show()
        else:
            plt.savefig('spherical_plot.png')


if __name__ == '__main__':
    argh.dispatch_command(main)
