#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import, unicode_literals


import os
from itertools import product

import numpy as np
import pytest

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL
from chemreac.integrate import run
from chemreac.serialization import load

"""
Tests the integration routine for the
chemical reaction system. (no diffusion)
"""

log = list(product([True, False], [True, False]))
@pytest.mark.parametrize("log", log)
def test_decay(log):
    # A -> B
    n = 2
    logy, logt = log
    k0 = 1.3
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], logy=logy, logt=logt)
    y0 = [3.0, 1.0]
    t0 = 5.0
    tend = 17.0
    nt = 42
    tout = np.linspace(t0, tend, nt+1)

    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(rd, y, t)
    if logy: yout = np.exp(yout)

    yref = np.array([y0[0]*np.exp(-k0*(tout-t0)),
                     y0[1]+y0[0]*(1-np.exp(-k0*(tout-t0)))]).transpose()
    import sys
    sys.stdout.write(yref-yout)
    assert np.allclose(yout, yref)


def test_ReactionDiffusion__bin_k_factor():
    # A -> B
    pass


@pytest.mark.parametrize("N", range(2,17))
def test_integrate__only_1_species_diffusion__mass_conservation(N):
    # Test that mass convervation is fulfilled wrt diffusion.
    x = np.linspace(0.1, 1.0, N+1)
    y0 = (x[0]/2+x[1:])**2

    geoms = (FLAT, SPHERICAL, CYLINDRICAL)

    for i, G in enumerate(geoms):
        sys = ReactionDiffusion(1, [], [], [], N=N, D=[0.02], x=x, geom=G)
        tout = np.linspace(0, 10.0, 50)
        yout, info = run(sys, y0, tout)
        if i == 0:
            yprim = yout
        elif i == 1:
            yprim = yout*(x[1:]**3-x[:-1]**3)
        else:
            yprim = yout*(x[1:]**2-x[:-1]**2)

        ybis = np.sum(yprim, axis=1)

        assert np.allclose(np.average(ybis), ybis)
