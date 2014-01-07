#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from chemreac.serialization import load

"""
Test chemical reaction system with 4 species.
(no diffusion)
"""

# A -> B               k1=0.05
# 2C + B -> D + B      k2=3.0

def _get_ref_f(sys, y0):
    k1, k2 = sys.k
    A, B, C, D = y0
    ref_f = np.array([-k1*A, k1*A, -2*k2*C*C*B, k2*C*C*B])
    return ref_f


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
    sys = load('four_species.json', N=N)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 0.0

    ref_f = _get_ref_f(sys, y0)
    fout = np.empty_like(ref_f)
    sys.f(t0, y0, fout)
    assert np.allclose(fout, ref_f)


def test_dense_jac_rmaj():
    N=1
    sys = load('four_species.json', N=N)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 0.0

    ref_J = _get_ref_J(sys, y0)
    Jout = np.empty_like(ref_J)
    sys.dense_jac_rmaj(t0, y0, Jout)
    assert np.allclose(Jout, ref_J)


def test_dense_jac_cmaj():
    N=1
    sys = load('four_species.json', N=N)
    y0 = np.array([1.3, 1e-4, 0.7, 1e-4])
    t0 = 0.0

    ref_J = _get_ref_J(sys, y0, order='F')
    Jout = np.empty_like(ref_J)
    sys.dense_jac_cmaj(t0, y0, Jout)
    assert np.allclose(Jout, ref_J)


if __name__ == '__main__':
    test_f()
    test_dense_jac_rmaj()
    test_dense_jac_cmaj()
