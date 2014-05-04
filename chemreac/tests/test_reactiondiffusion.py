#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from itertools import product
from math import exp

import numpy as np
import pytest

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL
from chemreac.util.banded import get_banded
from chemreac.util.grid import padded_centers, bounds, y_indices

np.set_printoptions(precision=3, linewidth=120)
TRUE_FALSE_PAIRS = list(product([True, False], [True, False]))


def test_autobinary():
    from chemreac.chemistry import (
        Reaction, ReactionSystem, mk_sn_dict_from_names
    )
    sbstncs = mk_sn_dict_from_names('AB')
    k = 3.0
    r1 = Reaction({'A': 2}, {'B': 1}, k=k)
    rsys = ReactionSystem([r1])
    rd = rsys.to_ReactionDiffusion(sbstncs)

    fout = np.empty((2,))
    rd.f(0.0, np.asarray([1.0, 37.0]), fout)
    assert np.allclose(fout, [-2*3.0, 3.0])


@pytest.mark.xfail
def test_ReactionDiffusion__f__wrong_fout_dimension():
    y0 = np.array([2.0, 3.0])
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k])
    fout = np.ones((1,))*99  # fout too small
    rd.f(0.0, y0, fout)


@pytest.mark.xfail
def test_ReactionDiffusion__too_few_species():
    # Ensure exception raised when referencing species indices > (n-1)
    # A -> B
    n = 1  # wrong: A, B makes 2
    with pytest.raises(RuntimeError):
        ReactionDiffusion(n, [[0]], [[1]], [k])


@pytest.mark.parametrize("N", [1, 3, 4])
def test_ReactionDiffusion__only_1_reaction(N):
    t0 = 3.0
    y0 = np.array([2.0, 3.0]*N)
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0])
    fout = np.ones((2*N,))*99
    rd.f(t0, y0, fout)

    for i in range(N):
        assert np.allclose(fout[i*2:(i+1)*2], np.array([-10.0, 10.0]))

    jout = np.zeros((2*N, 2*N))
    jref = np.zeros((2*N, 2*N))
    for i in range(N):
        jref[i*2,   i*2] = -k
        jref[i*2+1, i*2] = k
    rd.dense_jac_rmaj(t0, y0, jout)
    assert np.allclose(jout, jref)


def test_ReactionDiffusion__actv_1():
    y0 = np.array([2.0, 3.0, 7.0])
    k = 5.0
    # A + C -(+A)-> B + C
    rd = ReactionDiffusion(3, [[0, 0, 2]], [[1, 2]], [k], stoich_actv=[[0, 2]])
    fout = np.empty((3,))
    rd.f(0.0, y0, fout)
    r = k*y0[0]*y0[2]
    assert np.allclose(fout, [-2*r, r, 0])


def test_ReactionDiffusion__actv_2():
    y0 = np.array([2.0, 3.0, 9.0])
    k = 5.0
    # A + C --(+A+5*C)--> B
    rd = ReactionDiffusion(3, [[0, 0, 2, 2, 2, 2, 2, 2]], [[1]], [k],
                           stoich_actv=[[0, 2]])
    fout = np.empty((3,))
    rd.f(0.0, y0, fout)
    r = k*y0[0]*y0[2]
    assert np.allclose(fout, [-2*r, r, -6*r])


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_ReactionDiffusion__lrefl_3(log):
    # Diffusion without reaction
    # 3 bins
    t0 = 3.0
    logy, logt = log
    D = 17.0
    y0 = np.array([23.0, 27.0, 37.0])
    x = [5.0, 9.0, 13.0, 15.0]
    nstencil = 3
    xc = [3.0, 7.0, 11.0, 14.0, 16.0]

    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x, logy=logy,
                           nstencil=nstencil, logt=logt,
                           lrefl=True, rrefl=False)
    assert np.allclose(rd._xc, [3, 7, 11, 14, 16])
    fout = np.ones((3,))*99

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.f(t, y, fout)
    # In [7]: xlst=[0, 3, 7, 11, 14]

    # In [8]: print(finite_diff_weights(2, xlst[1:4], x0=xlst[2])[-1][-1])
    # [1/16, -1/8, 1/16]

    # In [9]: print(finite_diff_weights(2, xlst[2:5], x0=xlst[3])[-1][-1])
    # [1/14, -1/6, 2/21]

    # In [10]: print(finite_diff_weights(2, xlst[2:5], x0=xlst[4])[-1][-1])
    # [1/14, -1/6, 2/21]

    D_weight_ref = np.array([1/16, -1/8, 1/16, 1/14, -1/6,
                             2/21, 1/14, -1/6, 2/21])
    assert np.allclose(rd.D_weight, D_weight_ref)

    fref = np.array([
        1/16*y0[0] - 1/8*y0[0] + 1/16*y0[1],  # lrefl=True
        1/14*y0[0] - 1/6*y0[1] + 2/21*y0[2],
        1/14*y0[0] - 1/6*y0[1] + 2/21*y0[2],  # rrefl=False
    ])*D

    if logy:
        fref /= y0
    if logt:
        fref *= t0

    assert np.allclose(fout, fref)

    if logy:
        Jref = D*np.array([
            [
                -1/16*y0[1]/y0[0],
                1/16*y0[1]/y0[0],
                0
            ],
            [
                1/14*y0[0]/y0[1],
                -1/14*y0[0]/y0[1] - 2/21*y0[2]/y0[1],
                2/21*y0[2]/y0[1]
            ],
            [
                1/14*y0[0]/y0[2],
                -1/6*y0[1]/y0[2],
                -1/14*y0[0]/y0[2] + 1/6*y0[1]/y0[2]
            ]
        ])
    else:
        Jref = D*np.array([
            [1/16-1/8, 1/16, 0],
            [1/14, -1/6, 2/21],
            [1/14, -1/6, 2/21],
        ])

    if logt:
        Jref *= t0

    Jref[0,2] = 0  # Not perfect Jacobian
    Jref[2,0] = 0  # only nearest neighbour
    Jout = np.zeros_like(Jref)

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.dense_jac_rmaj(t, y, Jout)

    assert np.allclose(Jout, Jref)


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_ReactionDiffusion__rrefl_3(log):
    # Diffusion without reaction
    # 3 bins
    t0 = 3.0
    logy, logt = log
    D = 17.0
    y0 = np.array([23.0, 27.0, 37.0])
    x = [5.0, 9.0, 13.0, 15.0]
    nstencil = 3
    xc = [3.0, 7.0, 11.0, 14.0, 16.0]

    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x, logy=logy,
                           nstencil=nstencil, logt=logt,
                           lrefl=False, rrefl=True)
    assert np.allclose(rd._xc, [3, 7, 11, 14, 16])
    fout = np.ones((3,))*99

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.f(t, y, fout)
    # In [7]: xlst=[3, 7, 11, 14, 16]

    # In [8]: print(finite_diff_weights(2, xlst[1:4], x0=xlst[1])[-1][-1])
    # [1/14, -1/6, 2/21]

    # In [9]: print(finite_diff_weights(2, xlst[1:4], x0=xlst[2])[-1][-1])
    # [1/14, -1/6, 2/21]

    # In [10]: print(finite_diff_weights(2, xlst[2:5], x0=xlst[3])[-1][-1])
    # [2/15, -1/3, 1/5]

    D_weight_ref = np.array([1/14, -1/6, 2/21, 1/14, -1/6, 2/21, 2/15, -1/3, 1/5])
    assert np.allclose(rd.D_weight, D_weight_ref)

    fref = np.array([
        1/14*y0[0] - 1/6*y0[1] + 2/21*y0[2],  # lrefl=False
        1/14*y0[0] - 1/6*y0[1] + 2/21*y0[2],
        2/15*y0[1] - 1/3*y0[2] +  1/5*y0[2],  # rrefl=True
    ])*D

    if logy:
        fref /= y0
    if logt:
        fref *= t0

    assert np.allclose(fout, fref)

    if logy:
        Jref = D*np.array([
            [
                1/6*y0[1]/y0[0] - 2/21*y0[2]/y0[0],
                -1/6*y0[1]/y0[0],
                2/21*y0[2]/y0[0]
            ],
            [
                1/14*y0[0]/y0[1],
                -1/14*y0[0]/y0[1] - 2/21*y0[2]/y0[1],
                2/21*y0[2]/y0[1]
            ],
            [
                0,
                2/15*y0[1]/y0[2],
                -2/15*y0[1]/y0[2]
            ]
        ])
    else:
        Jref = D*np.array([
            [1/14, -1/6, 2/21],
            [1/14, -1/6, 2/21],
            [   0, 2/15, 1/5-1/3]
        ])

    if logt:
        Jref *= t0

    Jref[0,2] = 0  # Not perfect Jacobian
    Jref[2,0] = 0  # only nearest neighbour
    Jout = np.zeros_like(Jref)

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.dense_jac_rmaj(t, y, Jout)

    assert np.allclose(Jout, Jref)


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_ReactionDiffusion__lrefl_7(log):
    # Diffusion without reaction (7 bins)
    lrefl, rrefl = True, False
    N = 7
    t0 = 3.0
    logy, logt = log
    D = 17.0
    nstencil = 5
    nsidep = (nstencil-1)//2
    x = np.array([3, 5, 13, 17, 23, 25, 35, 37], dtype=np.float64)
    xc_ = padded_centers(x, nsidep)
    b = bounds(nstencil, N, lrefl=lrefl, rrefl=rrefl)

    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x, logy=logy,
                           logt=logt, N=N, nstencil=nstencil,
                           lrefl=lrefl, rrefl=rrefl)
    fout = np.ones((N,))*99

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.f(t, y, fout)
    le = nsidep if lrefl else 0
    D_weight_ref = np.array([
        finite_diff_weights(
            2, xc_[b[i][0]:b[i][1]], x0=xc_[le+i])[-1][-1]
        for i in range(N)])
    assert np.allclose(rd.D_weight, D_weight_ref.flatten())

    yi = y_indices(nstencil, N)

    fref = D*np.array([
        sum([rd.D_weight[i*nstencil+j]*y0[yi[j]]
             for j in range(b[i][0], b[i][1])])
        for i in range(N)
    ])

    if logy:
        fref /= y0
    if logt:
        fref *= t0

    assert np.allclose(fout, fref)

    if logy:
        Jref = D*np.array([
            [
                1/6*y0[1]/y0[0] - 2/21*y0[2]/y0[0],
                -1/6*y0[1]/y0[0],
                2/21*y0[2]/y0[0]
            ],
            [
                1/14*y0[0]/y0[1],
                -1/14*y0[0]/y0[1] - 2/21*y0[2]/y0[1],
                2/21*y0[2]/y0[1]
            ],
            [
                0,
                2/15*y0[1]/y0[2],
                -2/15*y0[1]/y0[2]
            ]
        ])
    else:
        def cb(i, j):
            if abs(i-j) > 1:
                return 0
            for k in range(nstencil):
                if yi[k] == j:
                    return rd.D_weight[i][k]
        Jref = np.fromfunction(cb, (N, N))
        Jref = D*np.array([
            [1/14, -1/6, 2/21],
            [1/14, -1/6, 2/21],
            [   0, 2/15, 1/5-1/3]
        ])

    if logt:
        Jref *= t0

    Jout = np.zeros_like(Jref)

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.dense_jac_rmaj(t, y, Jout)

    print('Jout')
    print(Jout)
    print('Jref')
    print(Jref)
    print('Jout-Jref')
    print(Jout-Jref)

    assert np.allclose(Jout, Jref)


@pytest.mark.parametrize("N", [1, 3, 4, 5])
def test_ReactionDiffusion__only_1_reaction__logy(N):
    # See <test_ReactionDiffusion__only_1_reaction__logy.png>
    t0 = 3.0
    y0 = np.array([2.0, 3.0]*N)
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0], logy=True)
    fout = np.ones((2*N,))*99
    rd.f(t0, np.log(y0), fout)

    for i in range(N):
        y0_ = y0[i*2:(i+1)*2]
        assert np.allclose(fout[i*2:(i+1)*2], [-k, k*y0_[0]/y0_[1]])

    jout = np.zeros((2*N, 2*N))
    jref = np.zeros((2*N, 2*N))
    for i in range(N):
        A = y0[i*2]
        B = y0[i*2+1]
        jref[i*2+1, i*2] = k/B*A
        jref[i*2+1, i*2+1] = -k/B*A
    rd.dense_jac_rmaj(t0, np.log(y0), jout)
    assert np.allclose(jout, jref)


@pytest.mark.parametrize("N", [1, 3, 4, 5])
def test_ReactionDiffusion__only_1_reaction__logy__logt(N):
    # See <test_ReactionDiffusion__only_1_reaction__logy_logt.png>
    t0 = 3.0
    y0 = np.array([2.0, 3.0]*N)
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0],
                           logy=True, logt=True)
    fout = np.ones((2*N,))*99
    rd.f(np.log(t0), np.log(y0), fout)

    for i in range(N):
        y0_ = y0[i*2:(i+1)*2]
        assert np.allclose(fout[i*2:(i+1)*2], [-k*t0, t0*k*y0_[0]/y0_[1]])


@pytest.mark.parametrize("N", [1, 3, 4, 5])
def test_ReactionDiffusion__only_1_reaction_bin_k_factor(N):
    y0 = np.concatenate([np.array([2.0, 3.0])/(x+1) for x in range(N)])
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0],
                           bin_k_factor=[[x + 1] for x in range(N)],
                           bin_k_factor_span=[1])
    fout = np.ones((2*N,))*99
    rd.f(0.0, y0, fout)

    for i in range(N):
        assert np.allclose(fout[i*2:(i+1)*2], np.array([-10.0, 10.0]))


@pytest.mark.parametrize("N", [1, 3, 4, 5])
def test_ReactionDiffusion__only_1_reaction_bin_k_factor_logy(N):
    y0 = np.concatenate([np.array([2.0, 3.0])/(x+1) for x in range(N)])
    k = 5.0
    # A -> B

    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0],
                           bin_k_factor=[[x+1] for x in range(N)],
                           bin_k_factor_span=[1], logy=True)
    fout = np.ones((2*N,))*99
    rd.f(0.0, np.log(y0), fout)

    def k_(bi):
        return k*(bi+1)

    for i in range(N):
        y0_ = y0[i*2:(i+1)*2]
        assert np.allclose(fout[i*2:(i+1)*2], [-k_(i), k_(i)*y0_[0]/y0_[1]])


@pytest.mark.parametrize("N", [1, 3, 4, 5])
def test_ReactionDiffusion__only_1_reaction_bin_k_factor_logy_logt(N):
    t0 = 3.0
    y0 = np.concatenate([np.array([2.0, 3.0])/(x+1) for x in range(N)])
    k = 5.0
    # A -> B

    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0],
                           bin_k_factor=[[x+1] for x in range(N)],
                           bin_k_factor_span=[1], logy=True, logt=True)
    fout = np.ones((2*N,))*99
    rd.f(np.log(t0), np.log(y0), fout)

    def k_(bi):
        return k*(bi+1)

    for i in range(N):
        y0_ = y0[i*2:(i+1)*2]
        assert np.allclose(fout[i*2:(i+1)*2],
                           [-k_(i)*t0, k_(i)*t0*y0_[0]/y0_[1]])


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_ReactionDiffusion__only_1_species_diffusion_3bins(log):
    # Diffusion without reaction
    # 3 bins
    t0 = 3.0
    logy, logt = log
    D = 17.0
    y0 = np.array([23.0, 27.0, 37.0])
    N = 3
    x = [5.0, 7.0, 13.0, 15.0]
    xc = [6.0, 10.0, 14.0]
    nstencil = 3
    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x, N=N, logy=logy,
                           logt=logt, lrefl=False, rrefl=False, nstencil=nstencil)
    assert np.allclose(rd._xc, [4, 6, 10, 14, 16])
    fout = np.ones((3,))*99

    w = [1/16, -1/8, 1/16]  # finite diff. weights for 2nd order deriv
    for i in range(N):
        assert np.allclose(rd.D_weight[i*nstencil:(i+1)*nstencil], w)
    J = D*(w[0]*y0[0] + w[1]*y0[1] + w[2]*y0[2])
    fref = np.array([J, J, J])

    if logy:
        fref /= y0
    if logt:
        fref *= t0

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.f(t, y, fout)

    assert np.allclose(fout, fref)

    # See <test_ReactionDiffusion__only_1_species_diffusion_2bins.png>
    jout = np.zeros((3, 3))
    if logy:
        jref = D*np.array([  # jref[i, k] = ...
            [w[k]*y0[k]/y0[i] if k != i else -1/y0[k]*sum([
                w[j]*y0[j] if j != k else 0 for j in range(3)
            ]) for k in range(3)] for i in range(3)
        ])
        jref[0, 2] = 0.0  # dense_jac_rmaj only computes banded approx.
        jref[2, 0] = 0.0  # same as above.
    else:
        jref = D*np.array([[w[k] if abs(k-i) < 2 else 0.0 for
                          k in range(3)] for i in range(3)])

    if logt:
        jref *= t0

    rd.dense_jac_rmaj(t, y, jout)

    assert np.allclose(jout, jref)

    jout_bnd = np.zeros((3, 3), order='F')
    rd.banded_packed_jac_cmaj(t, y, jout_bnd)
    jref_bnd = get_banded(jref, 1, 3)
    assert np.allclose(jout_bnd, jref_bnd)


def test_ReactionDiffusion__D_weight():
    x = np.array([2, 4, 6, 8, 10, 12, 14, 16], dtype=np.float64)
    # xc = x[:-1] + np.diff(x)/2
    rd = ReactionDiffusion(1, [], [], [], D=[1], x=x, nstencil=3)
    assert np.allclose(rd.D_weight, np.array([
        [1/4, -1/2, 1/4],
        [1/4, -1/2, 1/4],
        [1/4, -1/2, 1/4],
        [1/4, -1/2, 1/4],
        [1/4, -1/2, 1/4],
        [1/4, -1/2, 1/4],
        [1/4, -1/2, 1/4]
    ]).flatten())

    rd = ReactionDiffusion(1, [], [], [], D=[1], x=x, nstencil=5)
    assert np.allclose(rd.D_weight, np.array([
        [35/48, -13/6, 19/8, -7/6, 11/48],
        [11/48, -5/12, 1/8, 1/12, -1/48],
        [-1/48, 1/3, -5/8, 1/3, -1/48],
        [-1/48, 1/3, -5/8, 1/3, -1/48],
        [-1/48, 1/3, -5/8, 1/3, -1/48],
        [-1/48, 1/12, 1/8, -5/12, 11/48],
        [11/48, -7/6, 19/8, -13/6, 35/48]
    ]).flatten())


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_ReactionDiffusion__only_1_species_diffusion_7bins(log):
    # Diffusion without reaction
    N = 7
    nstencil = 5
    t0 = 3.0
    logy, logt = log
    D = 2.0
    y0 = np.array([12, 8, 11, 5, 7, 4, 9], dtype=np.float64)
    x = np.array([3, 5, 13, 17, 23, 25, 35, 37], dtype=np.float64)
    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x,
                           logy=logy, logt=logt, nstencil=nstencil,
                           lrefl=False, rrefl=False)
    weights = [
        [951/8800, -716/2475, 100/297, -75/352, 311/5400],
        [321/8800, -161/2475, 7/297, 3/352, -19/5400],
        [-39/8800, 109/2475, -127/1485, 87/1760, -19/5400],
        [-2/693, 38/675, -129/1100, 7/108, -1/1050],
        [0, 9/160, -7/72, 2/45, -1/288],
        [-8/1575, 9/400, 0, -19/450, 25/1008],
        [16/315, -9/32, 31/72, -13/45, 179/2016]
    ]
    assert np.allclose(rd.D_weight, np.array(weights).flatten())

    fout = np.ones((N,))*99
    fref = np.array([sum([D*weights[i][j]*y0[j+b[i][0]] for j
                          in range(nstencil)]) for i in range(N)])

    if logy:
        fref /= y0
    if logt:
        fref *= t0

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.f(t, y, fout)
    assert np.allclose(fout, fref)

    jref = np.zeros((N, N))
    jout = np.zeros((N, N))
    for i in range(N):
        for j in range(max(0, i-1), min(N, i+2)):
            if logy:
                if j == i+1 or j == i-1:
                    jref[i, j] = D*weights[i][j-b[i][0]]*y0[j]/y0[i]
            else:
                if i-1 <= j and j <= i+1:
                    jref[i, j] = D*weights[i][j-b[i][0]]

    if logt:
        jref *= t0

    rd.dense_jac_rmaj(t, y, jout)
    assert np.allclose(jout, jref)

    jout_bnd = np.zeros((3, N), order='F')
    rd.banded_packed_jac_cmaj(t, y, jout_bnd)
    jref_bnd = get_banded(jref, 1, N)
    assert np.allclose(jout_bnd, jref_bnd)


@pytest.mark.parametrize("geom", (FLAT, SPHERICAL, CYLINDRICAL))
def test_ReactionDiffusion__3_reactions_4_species_5_bins_k_factor(geom):
    # r[0]: A + B -> C
    # r[1]: D + C -> D + A + B
    # r[2]: B + B -> D
    pi = 3.141592653589793
    #              r[0]     r[1]       r[2]
    stoich_reac = [[0, 1], [2, 3],    [1, 1]]
    stoich_prod = [[2],    [0, 1, 3], [3]]
    stoich_actv = stoich_reac
    n = 4
    N = 5
    D = [2.5, 3.7, 5.11, 7.13]
    y0 = np.array([
        2.5, 1.2, 3.2, 4.3,
        2.7, 0.8, 1.6, 2.4,
        3.1, 0.3, 1.5, 1.8,
        3.3, 0.6, 1.6, 1.4,
        3.6, 0.9, 1.7, 1.2
    ])
    x = np.array([11.0, 13.0, 17.0, 23.0, 29.0, 37.0])
    k = [31.0, 37.0, 41.0]

    # (r[0], r[1]) modulations over bins
    bin_k_factor = [(i+3, i+4) for i in range(N)]
    bin_k_factor_span = [1, 1]

    rd = ReactionDiffusion(
        4, stoich_reac, stoich_prod, k, N, D, x,
        bin_k_factor=bin_k_factor,
        bin_k_factor_span=bin_k_factor_span, geom=geom)

    # Let's calculate f "by hand"
    if geom == FLAT:
        A = [1.0]*(N+1)
        Vincl = x
    elif geom == SPHERICAL:
        A = [4*pi*r**2 for r in x]
        Vincl = [4*pi*r**3/3 for r in x]
    elif geom == CYLINDRICAL:
        A = [2*pi*r for r in x]
        Vincl = [pi*r**2 for r in x]
    V = [Vincl[i+1]-Vincl[i] for i in range(N)]
    dx = np.diff((x[:-1]+x[1:])/2)

    def C_(si, bi):
        return y0[si+n*bi]

    def flux(si, bi):
        C = C_(si, bi)
        if bi > 0:  # previous
            Cp = y0[si+n*(bi-1)]
        else:
            Cp = C  # Neumann boundary conditions

        if bi < N-1:
            Cn = y0[si+n*(bi+1)]
        else:
            Cn = C  # Neumann boundary conditions

        f = 0.0
        if bi > 0:
            f -= D[si]*A[bi]*(C-Cp)/dx[bi-1]
        if bi < N-1:
            f += D[si]*A[bi+1]*(Cn-C)/dx[bi]
        return f

    r = [
        [k[0]*bin_k_factor[bi][0]*C_(0, bi)*C_(1, bi) for
         bi in range(N)],
        [k[1]*bin_k_factor[bi][1]*C_(3, bi)*C_(2, bi) for
         bi in range(N)],
        [k[2]*C_(1, bi)**2 for bi in range(N)],
    ]

    ref_f = np.array([[
        -r[0][bi] + r[1][bi] + flux(0, bi)/V[bi],
        -r[0][bi] + r[1][bi] - 2*r[2][bi] + flux(1, bi)/V[bi],
        r[0][bi] - r[1][bi] + flux(2, bi)/V[bi],
        r[2][bi] + flux(3, bi)/V[bi]
    ] for bi in range(N)]).flatten()

    # Compare to what is calculated using our C++ callback
    fout = np.empty(n*N)
    rd.f(0.0, y0, fout)
    assert np.allclose(fout, ref_f)

    # Now let's check that the Jacobian is correctly computed.
    def dfdC(bi, lri, lci):
        v = 0.0
        for ri in range(len(stoich_reac)):
            totl = (stoich_prod[ri].count(lri) -
                    stoich_reac[ri].count(lri))
            if totl == 0:
                continue
            actv = stoich_actv[ri].count(lci)
            if actv == 0:
                continue
            v += actv*totl*r[ri][bi]/C_(lci, bi)
        return v

    def jac_elem(ri, ci):
        bri, bci = ri // n, ci // n
        lri, lci = ri % n,  ci % n
        elem = 0.0
        if bri == bci:
            # on block diagonal
            elem += dfdC(bri, lri, lci)
            if lri == lci:
                if bri > 0:
                    elem -= D[lri]*A[bri]/dx[bri-1]/V[bri]
                if bri < N-1:
                    elem -= D[lri]*A[bri+1]/dx[bri]/V[bri]
        elif bri == bci - 1:
            if lri == lci:
                # super diagonal (right diffusion)
                elem += D[lri]*A[bri+1]/dx[bri]/V[bri]
        elif bri == bci + 1:
            if lri == lci:
                # sub diagonal (left diffusion)
                elem += D[lri]*A[bri]/dx[bri-1]/V[bri]
        return elem

    ref_j = np.zeros((n*N, n*N), order='C')
    for ri, ci in np.ndindex(n*N, n*N):
        ref_j[ri, ci] = jac_elem(ri, ci)

    jout_rmaj = np.zeros((n*N, n*N), order='C')
    rd.dense_jac_rmaj(0.0, y0, jout_rmaj)
    assert np.allclose(jout_rmaj, ref_j)

    jout_cmaj = np.zeros((n*N, n*N), order='F')
    rd.dense_jac_cmaj(0.0, y0, jout_cmaj)
    assert np.allclose(jout_cmaj, ref_j)

    ref_banded_j = get_banded(ref_j, n, N)

    jout_bnd_packed_cmaj = np.zeros((2*n+1, n*N), order='F')
    rd.banded_packed_jac_cmaj(0.0, y0, jout_bnd_packed_cmaj)

    plot = False
    if plot:
        import matplotlib.pyplot as plt
        from chemreac.util import coloured_spy
        fig = plt.figure()
        ax = fig.add_subplot(3, 1, 1)
        coloured_spy(ref_banded_j, ax=ax)
        plt.title('ref_banded_j')
        ax = fig.add_subplot(3, 1, 2)
        coloured_spy(jout_bnd_packed_cmaj, ax=ax)
        plt.title('jout_bnd_packed_cmaj')
        ax = fig.add_subplot(3, 1, 3)
        coloured_spy(ref_banded_j-jout_bnd_packed_cmaj, ax=ax)
        plt.title('diff')
        plt.show()

    assert np.allclose(jout_bnd_packed_cmaj, ref_banded_j)

    jout_bnd_padded_cmaj = np.zeros((3*n+1, n*N), order='F')
    rd.banded_padded_jac_cmaj(0.0, y0, jout_bnd_padded_cmaj)
    assert np.allclose(jout_bnd_padded_cmaj[n:, :], ref_banded_j)
