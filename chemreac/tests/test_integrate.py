#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import, unicode_literals


import os
from itertools import product

import pytest
import numpy as np

from chemreac.serialization import load
from chemreac.integrate import run

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL


"""
Tests the integration routine for the
chemical reaction system. (no diffusion)
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0

json_path, blessed_path = map(
    lambda x: os.path.join(os.path.dirname(__file__), x),
    ['four_species.json', 'four_species_blessed.txt']
)


combos = list(product([True, False], [True, False], range(1,4), [FLAT, SPHERICAL, CYLINDRICAL]))
@pytest.mark.parametrize("combo", combos)
def test_integrate(combo):
    logy, logt, N, geom = combo
    sys = load(json_path, N=N, logy=logy, logt=logt, geom=geom)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4]*N)

    ref = np.genfromtxt(blessed_path)
    ref_t = ref[:,0]
    ref_y = ref[:,1:5]

    t0 = 3.0
    tend=10.0+t0
    nt=100
    tout = np.linspace(t0, tend, nt+1)

    y = np.log(y0) if logy else y0
    if logt:
        tout = np.log(tout)
    tout, yout, info = run(sys, y, tout)
    if logy: yout = np.exp(yout)
    if logt: tout = np.exp(tout)
    tout -= t0

    assert np.allclose(tout, ref_t)
    assert np.allclose(yout[:,:4], ref_y, atol=1e-5)


@pytest.mark.parametrize("N", range(2,17))
def test_integrate__only_1_species_diffusion__mass_conservation(N):
    # Test that mass convervation is fulfilled wrt diffusion.
    x = np.linspace(0.1, 1.0, N+1)
    y0 = (x[0]/2+x[1:])**2

    geoms = (FLAT, SPHERICAL, CYLINDRICAL)

    for i, G in enumerate(geoms):
        sys = ReactionDiffusion(1, [], [], [], N=N, D=[0.02], x=x, geom=G)
        tout = np.linspace(0, 10.0, 50)
        tout, yout, info = run(sys, y0, tout)
        if i == 0:
            yprim = yout
        elif i == 1:
            yprim = yout*(x[1:]**3-x[:-1]**3)
        else:
            yprim = yout*(x[1:]**2-x[:-1]**2)

        ybis = np.sum(yprim, axis=1)

        assert np.allclose(np.average(ybis), ybis)
