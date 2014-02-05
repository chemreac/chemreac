#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from itertools import product

import numpy as np
import pytest

from chemreac.serialization import load

"""
Test chemical reaction system with 4 species.
(no diffusion)

tests:
chemreac.serialization.load
chemreac.PyReactionDiffusion.f
chemreac.PyReactionDiffusion.dense_jac_rmaj
chemreac.PyReactionDiffusion.dense_jac_cmaj

See:
<four_species_f_jac.png>
<four_species_f_jac_logy.png>
<four_species_f_jac_logt.png>
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0

JSON_PATH = os.path.join(os.path.dirname(__file__), 'four_species.json')
TRUE_FALSE_PAIRS = list(product([True, False], [True, False]))


def _get_ref_f(sys, t0, y0, logy, logt):
    k1, k2 = sys.k
    A, B, C, D = y0

    ref_f = np.array([-k1*A, k1*A, -2*k2*C*C*B, k2*C*C*B])
    if logy:
        ref_f /= y0
    if logt:
        ref_f *= t0
    return ref_f


def _get_ref_J(sys, t0, y0, logy, logt, order='C'):
    k1, k2 = sys.k
    A, B, C, D = y0
    if logy:
        f1,f2,f3,f4 = _get_ref_f(sys, t0, y0, logy, logt)
        ref_J = np.array(
            [[   0,   0,    0,   0],
             [  f2, -f2,    0,   0],
             [   0,  f3,   f3,   0],
             [   0,  f4, 2*f4, -f4]],
            order=order)
    else:
        ref_J = np.array(
            [[-k1,         0,           0, 0],
             [ k1,         0,           0, 0],
             [  0, -2*k2*C*C, -2*2*k2*C*B, 0],
             [  0,    k2*C*C,    2*k2*C*B, 0]],
            order=order)
        if logt:
            ref_J *= t0

    return ref_J


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_f(log):
    logy, logt = log
    N=1
    sys = load(JSON_PATH, N=N, logy=logy, logt=logt)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 42.0

    ref_f = _get_ref_f(sys, t0, y0, logy, logt)
    fout = np.empty_like(ref_f)

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    sys.f(t, y, fout)
    assert np.allclose(fout, ref_f)


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_dense_jac_rmaj(log):
    logy, logt = log
    N=1
    sys = load(JSON_PATH, N=N, logy=logy, logt=logt)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 42.0

    ref_J = _get_ref_J(sys, t0, y0, logy, logt)
    Jout = np.empty_like(ref_J)

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    sys.dense_jac_rmaj(t, y, Jout)

    assert np.allclose(Jout, ref_J)


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_dense_jac_cmaj(log):
    logy, logt = log
    N=1
    sys = load(JSON_PATH, N=N, logy=logy, logt=logt)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 42.0

    ref_J = _get_ref_J(sys, t0, y0, logy, logt, order='F')
    Jout = np.empty_like(ref_J)

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    sys.dense_jac_cmaj(t, y, Jout)

    assert np.allclose(Jout, ref_J)
