#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argh
import numpy as np

from reactiondiffusion.methods import (
    demo_methods, std_methods, EulerBackward,
    RefScipy, EulerBackward_dgbsv, EulerBackward_dgesv
)
from reactiondiffusion.serialization import load

"""
Demo of chemical reaction diffusion system.
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0

def main(nt=50, tend=10.0, use_extrapol_jac=False, stepper='DormandPrince45', implementation='py', N=None, demo=True):
    if implementation == 'py':
        from reactiondiffusion.chem import ReactionDiffusion
    elif implementation == 'cy':
        from reactiondiffusion.cy_chem import ReactionDiffusion
    elif implementation == 'cpp':
        from reactiondiffusion.cpp_chem_wrapper import PyReactionDiffusion as ReactionDiffusion
    else:
        raise NotImplementedError

    sys = load('four_species.json', RD=ReactionDiffusion, N=N)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4]+\
                  [0.1, 0.1, 0.1, 0.1]*(sys.N-1))


    t0 = 0.0


    if demo:
        method = filter(
            lambda x: x.__name__.lower() == stepper.lower(),
            std_methods)[0]

        demo_methods([method, RefScipy],
                     #std_methods, #[EulerBackward_dgbsv, RefScipy],
                     sys, y0, t0, tend, nt,
                     show=True, use_extrapol_jac=use_extrapol_jac)
    else:
        # Used for testing
        assert sys.N == 1 # Leave diffusion of of the picture for now
        k1, k2 = sys.k
        A, B, C, D = y0
        ref_f = np.array([-k1*A, k1*A, -2*k2*C*C*B, k2*C*C*B])

        ref_J = np.array([[-k1,         0,         0, 0],
                          [ k1,         0,         0, 0],
                          [  0, -2*k2*C*C, -2*2*k2*C*B, 0],
                          [  0,    k2*C*C,   2*k2*C*B, 0]])
        ref_J
        return sys, y0, t0, ref_f, ref_J

if __name__ == '__main__':
    argh.dispatch_command(main)
