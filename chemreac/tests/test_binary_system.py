#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division, print_function)

import os
from itertools import product

import numpy as np
import pytest

from chemreac import ReactionDiffusion
from chemreac.integrate import run
from chemreac.serialization import load
from chemreac.chemistry import mk_sn_dict_from_names, Reaction, ReactionSystem
from chemreac.util.testing import slow, veryslow

from .test_reactiondiffusion import _test_f_and_dense_jac_rmaj

"""
Test chemical reaction system with 4 species.
(no diffusion)

A + B -> C           k1=0.3

tests:
* chemreac.serialization.load
* chemreac.PyReactionDiffusion.f
* chemreac.PyReactionDiffusion.dense_jac_rmaj
* chemreac.PyReactionDiffusion.dense_jac_cmaj
* chemreac.integrate.run
* chemreac.chemistry.Reaction
* chemreac.chemistry.mk_sn_dict_from_names
* chemreac.chemistry.ReactionSystem
"""

TR_FLS = (True, False)

JSON_PATH, BLESSED_PATH = map(
    lambda x: os.path.join(os.path.dirname(__file__), x),
    ['binary_system.json', 'binary_system_blessed.txt']
)
TRUE_FALSE_PAIRS = list(product(TR_FLS, TR_FLS))
TRUE_FALSE_TRIPLES = list(product(TR_FLS, TR_FLS, TR_FLS))
C0 = [5.0, 11.0, 7.0]
D = [0.1, 0.2, 0.3]


def test_serialization():
    rd = load(JSON_PATH)
    assert rd.n == 3
    assert rd.stoich_active == [[0, 1]]
    assert rd.stoich_prod == [[2]]
    assert list(rd.k) == [0.3]
    assert rd.D.tolist() == [0.1, 0.2, 0.3]


def _get_ref_f(rd, t0, y0, logy, logt, use_log2=False):
    k = rd.k[0]
    A, B, C = y0

    ref_f = np.array([-k*A*B, -k*A*B, k*A*B])
    if logy:
        ref_f /= y0
        if not logt and use_log2:
            ref_f /= np.log(2)
    if logt:
        ref_f *= t0
        if not logy and use_log2:
            ref_f *= np.log(2)
    return ref_f


def _get_ref_J(rd, t0, y0, logy, logt, order='C', use_log2=False):
    k = rd.k[0]
    A, B, C = y0
    if logy:
        f1, f2, f3 = _get_ref_f(rd, t0, y0, logy, logt)
        if logt:
            A *= t0
            B *= t0
        ref_J = np.array(
            [[0,    -k*B,   0],
             [-k*A,    0,   0],
             [f3,     f3, -f3]],
            order=order)
    else:
        ref_J = np.array(
            [[-k*B, -k*A, 0],
             [-k*B, -k*A, 0],
             [k*B,   k*A, 0]],
            order=order)
        if logt:
            ref_J *= t0
    if logt and use_log2:
        ref_J *= np.log(2)

    return ref_J


@pytest.mark.parametrize("log", TRUE_FALSE_TRIPLES[::-1])
def test_f(log):
    logy, logt, use_log2 = log
    N = 1
    rd = load(JSON_PATH, N=N, logy=logy, logt=logt, use_log2=use_log2)
    y0 = np.array(C0)
    t0 = 42.0
    ref_f = _get_ref_f(rd, t0, y0, logy, logt, use_log2=use_log2)
    fout = rd.alloc_fout()

    y = rd.logb(y0) if logy else y0
    t = rd.logb(t0) if logt else t0
    rd.f(t, y, fout)
    assert np.allclose(fout, ref_f)


@pytest.mark.parametrize("log", TRUE_FALSE_TRIPLES)
def test_dense_jac_rmaj(log):
    logy, logt, use_log2 = log
    N = 1
    rd = load(JSON_PATH, N=N, logy=logy, logt=logt, use_log2=use_log2)
    y0 = np.array(C0)
    t0 = 42.0

    ref_J = _get_ref_J(rd, t0, y0, logy, logt, use_log2=use_log2)
    Jout = np.empty_like(ref_J)

    y = rd.logb(y0) if logy else y0
    t = rd.logb(t0) if logt else t0
    rd.dense_jac_rmaj(t, y, Jout)

    assert np.allclose(Jout, ref_J)


@pytest.mark.parametrize("log", TRUE_FALSE_TRIPLES)
def test_dense_jac_cmaj(log):
    logy, logt, use_log2 = log
    N = 1
    rd = load(JSON_PATH, N=N, logy=logy, logt=logt, use_log2=use_log2)
    y0 = np.array(C0)
    t0 = 42.0

    ref_J = _get_ref_J(rd, t0, y0, logy, logt, order='F', use_log2=use_log2)
    Jout = np.empty_like(ref_J)

    y = rd.logb(y0) if logy else y0
    t = rd.logb(t0) if logt else t0
    rd.dense_jac_cmaj(t, y, Jout)

    assert np.allclose(Jout, ref_J)


def test_chemistry():
    sbstncs = mk_sn_dict_from_names('ABC', D=[.1, .2, .3])
    r1 = Reaction({'A': 1, 'B': 1}, {'C': 1}, 0.3)
    rsys = ReactionSystem([r1], sbstncs)
    rd = ReactionDiffusion.from_ReactionSystem(rsys)
    serialized_rd = load(JSON_PATH)
    assert rd.stoich_active == serialized_rd.stoich_active
    assert rd.stoich_prod == serialized_rd.stoich_prod
    assert rd.stoich_inact == serialized_rd.stoich_inact
    assert np.allclose(rd.k, serialized_rd.k)
    assert np.allclose(rd.D, serialized_rd.D)


COMBOS = list(product(TR_FLS, TR_FLS, [1], 'f', TR_FLS))
SLOW_COMBOS = list(product(TR_FLS, TR_FLS, [3], 'fcs', TR_FLS))
VERYSLOW_COMBOS = list(product(TR_FLS, TR_FLS, [4, 5], 'fcs', TR_FLS))


def _test_integrate(params):
    logy, logt, N, geom, use_log2 = params
    rd = load(JSON_PATH, N=N, logy=logy, logt=logt, geom=geom, use_log2=use_log2)
    assert rd.geom == geom
    y0 = np.array(C0*N)

    ref = np.genfromtxt(BLESSED_PATH)
    ref_t = ref[:, 0]
    ref_y = ref[:, 1:rd.n+1]

    t0 = 3.0
    tend = 5.0+t0
    nt = 137
    tout = np.linspace(t0, tend, nt+1)
    assert np.allclose(tout-t0, ref_t)

    _test_f_and_dense_jac_rmaj(rd, t0, rd.logb(y0) if logy else y0)
    integr = run(rd, y0, tout, with_jacobian=True)
    for i in range(N):
        assert np.allclose(integr.Cout[:, i, :], ref_y, atol=1e-4)


@pytest.mark.parametrize("params", COMBOS)
def test_integrate(params):
    _test_integrate(params)


@slow
@pytest.mark.parametrize("params", SLOW_COMBOS)
def test_integrate_slow(params):
    _test_integrate(params)


@veryslow
@pytest.mark.parametrize("params", VERYSLOW_COMBOS)
def test_integrate_veryslow(params):
    _test_integrate(params)
