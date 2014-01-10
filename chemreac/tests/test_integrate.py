#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import pytest
import numpy as np

from chemreac.serialization import load
from chemreac.integrate import run

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
@pytest.mark.parametrize("N", range(1,5))
def test_integrate(N):
    sys = load(json_path, N=N)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4]*N)

    ref = np.genfromtxt(blessed_path)
    ref_t = ref[:,0]
    ref_y = ref[:,1:5]

    t0 = 0.0
    tend=10.0
    nt=100
    tout, yout, info = run(sys, y0, t0, tend, nt)

    assert np.allclose(tout, ref_t)
    assert np.allclose(yout[:,:4], ref_y, atol=1e-3)
