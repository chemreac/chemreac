#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from itertools import product

import numpy as np
import pytest

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL
from chemreac.integrate import run
from chemreac.serialization import load
from chemreac.chemistry import mk_sn_dict_from_names, Reaction, ReactionSystem

from test_reactiondiffusion import _test_dense_jac_rmaj

"""
Test chemical reaction system with 4 species.
(no diffusion)

A -> B               k1=0.05
2C + B -> D + B      k2=3.0

tests:
* chemreac.serialization.load
* chemreac.PyReactionDiffusion.f
* chemreac.PyReactionDiffusion.dense_jac_rmaj
* chemreac.PyReactionDiffusion.dense_jac_cmaj
* chemreac.integrate.run
* chemreac.chemistry.Reaction
* chemreac.chemistry.mk_sn_dict_from_names
* chemreac.chemistry.ReactionSystem

See:
<four_species_f_jac.png>
<four_species_f_jac_logy.png>
<four_species_f_jac_logt.png>
"""

np.set_printoptions(precision=3, linewidth=180)

JSON_PATH, BLESSED_PATH = map(
    lambda x: os.path.join(os.path.dirname(__file__), x),
    ['four_species.json', 'four_species_blessed.txt']
)
TRUE_FALSE_PAIRS = list(product([True, False], [True, False]))

Ns = [1] # change to [1, 3, 7]
combos = list(product([True, False], [True, False], [1], [FLAT, SPHERICAL, CYLINDRICAL]))
@pytest.mark.parametrize("combo", combos)
def test_integrate(combo):
    logy, logt, N, geom = combo
    sys = load(JSON_PATH, N=N, logy=logy, logt=logt, geom=geom)

    y0 = np.array([1.3, 1e-4, 0.7, 1e-4]*N)

    ref = np.genfromtxt(BLESSED_PATH)
    ref_t = ref[:,0]
    ref_y = ref[:,1:5]

    t0 = 3.0
    tend=10.0+t0
    nt=100
    tout = np.linspace(t0, tend, nt+1)
    assert np.allclose(tout-t0, ref_t)

    y = np.log(y0) if logy else y0
    if logt:
        tout = np.log(tout)
    yout, info = run(sys, y, tout)
    if logy: yout = np.exp(yout)

    for i in range(N):
        assert np.allclose(yout[:, i, :], ref_y, atol=1e-5)


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

    print(Jout)
    print(ref_J)
    print(Jout-ref_J)
    assert np.allclose(Jout, ref_J)


def test_chemistry():
    sbstncs = mk_sn_dict_from_names('ABCD', D=[0.1, 0.2, 0.3, 0.4])
    r1 = Reaction({'A': 1}, {'B': 1}, k=0.05)
    r2 = Reaction({'B': 1, 'C': 2}, {'D': 1, 'B': 1}, k=3.0)
    rsys = ReactionSystem([r1, r2])
    rd = rsys.to_ReactionDiffusion(sbstncs)
    # how to compare equality of say: json loaded sys?
    # specie indices can be permuted in 4*3*2*1 = 24 ways
    # ...solution: canonical representation is alphabetically sorted on
    #              Substance.name
    serialized_rd = load(JSON_PATH)
    assert rd.stoich_reac == serialized_rd.stoich_reac
    assert rd.stoich_prod == serialized_rd.stoich_prod
    assert rd.stoich_actv == serialized_rd.stoich_actv
    assert np.allclose(rd.k, serialized_rd.k)
    assert np.allclose(rd.D, serialized_rd.D)

def test_multi_compartment():
    rsys = load(JSON_PATH, N=3, lrefl=True, rrefl=True)
    _test_dense_jac_rmaj(rsys, 1.0, np.asarray([1.3, 1e-4, 0.7, 1e-4]*rsys.N).flatten())


if __name__ == '__main__':
    # Some old debug code which shows how to inspect using LaTeX:
    from chemreac.symbolic import SymRD
    import sympy as sp
    rd = load(JSON_PATH, N=3, lrefl=True, rrefl=True)
    print(rd.x)
    print(rd.stoich_reac)
    print(rd.stoich_prod)
    srd = SymRD.from_rd(rd)
    print('===   f ====')
    print('\n'.join(map(sp.latex, srd._f)))
    print('=== dfdy ====')
    for ri, row in enumerate(srd.jacobian.tolist()):
        for ci, expr in enumerate(row):
            if expr == 0:
                continue
            print(ri, ci, sp.latex(expr))
