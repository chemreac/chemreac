#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np

from chemreac.serialization import load

"""
Test chemical reaction system with 4 species.
(no diffusion)

tests:
chemreac.serialization.load
chemreac.PyReactionDiffusion.f
chemreac.PyReactionDiffusion.dense_jac_rmaj
chemreac.PyReactionDiffusion.dense_jac_cmaj

"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0

json_path = os.path.join(os.path.dirname(__file__), 'four_species.json')

def _get_ref_J(sys, y0, order='C'):
    k1, k2 = sys.k
    A, B, C, D = y0
    ref_J = np.array([[-k1,         0,           0, 0],
                      [ k1,         0,           0, 0],
                      [  0, -2*k2*C*C, -2*2*k2*C*B, 0],
                      [  0,    k2*C*C,    2*k2*C*B, 0]],
                     order=order)
    return ref_J


def test_f():
    N=1
    sys = load(json_path, N=N)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 0.0

    k1, k2 = sys.k
    A, B, C, D = y0
    ref_f = np.array([-k1*A, k1*A, -2*k2*C*C*B, k2*C*C*B])
    fout = np.empty_like(ref_f)
    sys.f(t0, y0, fout)
    assert np.allclose(fout, ref_f)


def test_dense_jac_rmaj():
    N=1
    sys = load(json_path, N=N)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 0.0

    ref_J = _get_ref_J(sys, y0)
    Jout = np.empty_like(ref_J)
    sys.dense_jac_rmaj(t0, y0, Jout)
    assert np.allclose(Jout, ref_J)


def test_dense_jac_cmaj():
    N=1
    sys = load(json_path, N=N)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 0.0

    ref_J = _get_ref_J(sys, y0, order='F')
    Jout = np.empty_like(ref_J)
    sys.dense_jac_cmaj(t0, y0, Jout)
    assert np.allclose(Jout, ref_J)
