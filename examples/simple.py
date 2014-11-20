#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

import argh
import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import Integration


def main(logy=False, logt=False):
    # A -> B
    n = 2
    k0 = 0.13
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], logy=logy, logt=logt)
    y0 = [3.0, 1.0]
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt+1)

    Cref = np.array([y0[0]*np.exp(-k0*(tout-t0)),
                     y0[1] + y0[0]*(1 - np.exp(-k0*(tout-t0)))]).transpose()

    # scipy
    integr1 = Integration('scipy', rd, np.asarray(y0), np.asarray(tout))
    assert np.allclose(integr1.Cout[:, 0, :], Cref)

    # sundials
    integr2 = Integration('cvode_direct', rd, np.asarray(y0), np.asarray(tout),
                          atol=[1e-8, 1e-8], rtol=1e-8, method='bdf')
    assert np.allclose(integr2.Cout[:, 0, :], Cref)


if __name__ == '__main__':
    argh.dispatch_command(main)
