#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

from itertools import product
from math import exp

import numpy as np
import pytest

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL

np.set_printoptions(precision=3, linewidth=120)
TRUE_FALSE_PAIRS = list(product([True, False], [True, False]))


def test_autobinary():
    from chemreac.chemistry import Reaction, ReactionSystem, mk_sn_dict_from_names
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
    fout = np.ones((1,))*99 # fout too small
    rd.f(0.0, y0, fout)


@pytest.mark.xfail
def test_ReactionDiffusion__too_few_species():
    # Ensure exception raised when referencing species indices > (n-1)
    y0 = np.array([2.0, 3.0])
    k = 5.0
    # A -> B
    n = 1 # wrong: A, B makes 2
    rd = ReactionDiffusion(1, [[0]], [[1]], [k])


@pytest.mark.parametrize("N", [1,3,4])
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

    jout = np.zeros((2*N,2*N))
    jref = np.zeros((2*N,2*N))
    for i in range(N):
        jref[i*2,  i*2] = -k
        jref[i*2+1,i*2] =  k
    rd.dense_jac_rmaj(t0, y0, jout)
    assert np.allclose(jout, jref)


def test_ReactionDiffusion__actv_1():
    y0 = np.array([2.0, 3.0, 7.0])
    k = 5.0
    # A + C -(+A)-> B + C
    rd = ReactionDiffusion(3, [[0,0,2]], [[1,2]], [k], stoich_actv=[[0,2]])
    fout = np.empty((3,))
    rd.f(0.0, y0, fout)
    r = k*y0[0]*y0[2]
    assert np.allclose(fout, [-2*r, r, 0])


def test_ReactionDiffusion__actv_2():
    y0 = np.array([2.0, 3.0, 9.0])
    k = 5.0
    # A + C --(+A+5*C)--> B
    rd = ReactionDiffusion(3, [[0,0,2,2,2,2,2,2]], [[1]], [k], stoich_actv=[[0,2]])
    fout = np.empty((3,))
    rd.f(0.0, y0, fout)
    r = k*y0[0]*y0[2]
    assert np.allclose(fout, [-2*r, r, -6*r])


@pytest.mark.parametrize("N", [1,3,4,5])
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

    jout = np.zeros((2*N,2*N))
    jref = np.zeros((2*N,2*N))
    for i in range(N):
        A = y0[i*2]
        B = y0[i*2+1]
        jref[i*2+1,i*2] = k/B*A
        jref[i*2+1,i*2+1] = -k/B*A
    rd.dense_jac_rmaj(t0, np.log(y0), jout)
    assert np.allclose(jout, jref)


@pytest.mark.parametrize("N", [1,3,4,5])
def test_ReactionDiffusion__only_1_reaction__logy__logt(N):
    # See <test_ReactionDiffusion__only_1_reaction__logy_logt.png>
    t0 = 3.0
    y0 = np.array([2.0, 3.0]*N)
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0], logy=True, logt=True)
    fout = np.ones((2*N,))*99
    rd.f(np.log(t0), np.log(y0), fout)

    for i in range(N):
        y0_ = y0[i*2:(i+1)*2]
        assert np.allclose(fout[i*2:(i+1)*2], [-k*t0, t0*k*y0_[0]/y0_[1]])


@pytest.mark.parametrize("N", [1,3,4,5])
def test_ReactionDiffusion__only_1_reaction_bin_k_factor(N):
    y0 = np.concatenate([np.array([2.0, 3.0])/(x+1) for x in range(N)])
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0],
                           bin_k_factor = [[x+1] for x in range(N)],
                           bin_k_factor_span=[1])
    fout = np.ones((2*N,))*99
    rd.f(0.0, y0, fout)

    for i in range(N):
        assert np.allclose(fout[i*2:(i+1)*2], np.array([-10.0, 10.0]))


@pytest.mark.parametrize("N", [1,3,4,5])
def test_ReactionDiffusion__only_1_reaction_bin_k_factor_logy(N):
    y0 = np.concatenate([np.array([2.0, 3.0])/(x+1) for x in range(N)])
    k = 5.0
    # A -> B

    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0],
                           bin_k_factor = [[x+1] for x in range(N)],
                           bin_k_factor_span=[1], logy=True)
    fout = np.ones((2*N,))*99
    rd.f(0.0, np.log(y0), fout)

    def k_(bi):
        return k*(bi+1)

    for i in range(N):
        y0_ = y0[i*2:(i+1)*2]
        assert np.allclose(fout[i*2:(i+1)*2], [-k_(i), k_(i)*y0_[0]/y0_[1]])


@pytest.mark.parametrize("N", [1,3,4,5])
def test_ReactionDiffusion__only_1_reaction_bin_k_factor_logy_logt(N):
    t0 = 3.0
    y0 = np.concatenate([np.array([2.0, 3.0])/(x+1) for x in range(N)])
    k = 5.0
    # A -> B

    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0],
                           bin_k_factor = [[x+1] for x in range(N)],
                           bin_k_factor_span=[1], logy=True, logt=True)
    fout = np.ones((2*N,))*99
    rd.f(np.log(t0), np.log(y0), fout)

    def k_(bi):
        return k*(bi+1)

    for i in range(N):
        y0_ = y0[i*2:(i+1)*2]
        assert np.allclose(fout[i*2:(i+1)*2], [-k_(i)*t0, k_(i)*t0*y0_[0]/y0_[1]])



def _get_banded(A, n, N):
    B = np.zeros((2*n+1, n*N))
    for ri in range(n*N):
        for ci in range(max(0, ri-n), min(n*N, ri+n+1)):
            B[n+ri-ci, ci] = A[ri, ci]
    return B

def test__get_banded():
    A = np.array([[2.0, 3.0],
                  [5.0, 7.0]])
    B = _get_banded(A, 1, 2)
    B_ref = np.array([
        [0.0, 3.0],
        [2.0, 7.0],
        [5.0, 0.0]
    ])
    assert np.allclose(B, B_ref)


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_ReactionDiffusion__only_1_species_diffusion_3bins(log):
    # Diffusion without reaction
    # 3 bins
    t0 = 3.0
    logy, logt = log
    D = 17.0
    y0 = np.array([23.0, 27.0, 37.0])
    y0 = np.array([1.0, 2.0, 1.0]) ## DEBUG
    x = [5.0, 7.0, 13.0, 15.0]
    xc = [6.0, 10.0, 14.0]
    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x, logy=logy, logt=logt)
    fout = np.ones((3,))*99

    w = [1/16, -1/8, 1/16]  # finite diff. weights for 2nd order deriv
    J = -D*(w[0]*y0[0] + w[1]*y0[1] + w[2]*y0[2])
    fref = [J, J, J]

    if logy:
        fref /= y0
    if logt:
        fref *= t0

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.f(t, y, fout)


    print(w)
    print(rd.D_weight) # DEBUG
    print(fout) # DEBUG
    print(fref) # DEBUG
    assert np.allclose(fout, fref)


    # See <test_ReactionDiffusion__only_1_species_diffusion_2bins.png>
    jout = np.zeros((3,3))
    if logy:
        jref = np.array([
            [-D*w[k]*exp(y0[k]-y0[i]) if k !=i else exp(-y0[k])*D*sum(
                [w[j]*exp(y0[j]-y0[k]) if j != k else 0 for j in range(3)])
             for k in range(3)] for i in range(3)
        ])
    else:
        jref = np.array([
            [-D*w[k]*y0[k] for k in range(3)] for i in range(3)
        ])

    if logt:
        jref *= t0

    rd.dense_jac_rmaj(t, y, jout)

    assert np.allclose(jout, jref)

    jout_bnd = np.zeros((3,3), order='F')
    rd.banded_packed_jac_cmaj(t, y, jout_bnd)
    jref_bnd = _get_banded(jref, 1, 3)
    assert np.allclose(jout_bnd, jref_bnd)


@pytest.mark.parametrize("log", TRUE_FALSE_PAIRS)
def test_ReactionDiffusion__only_1_species_diffusion_99bins(log):
    # Diffusion without reaction
    # 3 bins
    # See
    # <test_ReactionDiffusion__only_1_species_diffusion_3bins.png>
    # <test_ReactionDiffusion__only_1_species_diffusion_3bins_logy.png>
    # <only_1_species_diffusion_3bins_logy_formulae.png>
    N = 3
    t0 = 3.0
    logy, logt = log
    D = 2.0
    y0 = np.array([12., 8., 11.])
    x = [3., 5., 13., 17.]
    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x, logy=logy, logt=logt)
    fout = np.ones((N,))*99
    V = [x[i+1]-x[i] for i in range(N)]
    l = np.diff(x)
    dx = np.array([0.5*(x[i+2]-x[i]) for i in range(N-1)])
    dC = np.array([y0[i+1]-y0[i] for i in range(N-1)])
    J = -D*dC/dx
    J_ = np.pad(J, ((1,1),), mode='constant')
    fref = np.array([(J_[i]-J_[i+1])/l[i] for i in range(N)])
    assert np.allclose(fref, [-4./5, 13./40, -1./4])
    if logy:
        fref /= y0
        a,b,c = np.log(y0)
        assert np.allclose(fref, np.array([
            D/dx[0]/l[0]*(exp(b-a)-1),
            D/l[1]*((exp(c-b)-1)/dx[1]-(1-exp(a-b))/dx[0]),
            D/l[2]/dx[1]*(exp(b-c)-1)
         ]))
    if logt:
        fref *= t0

    y = np.log(y0) if logy else y0
    t = np.log(t0) if logt else t0
    rd.f(t, y, fout)
    assert np.allclose(fout, fref)

    jout = np.zeros((N,N))
    if logy:
        unlogt = t0 if logt else 1.0
        jref = np.array([
            [-fref[0]/unlogt-D/dx[0]/V[0],
             D/dx[0]/V[0]*exp(b-a),
             0],
            [D/dx[0]/V[1]*exp(a-b),
             -fref[1]/unlogt-D/V[1]*(1/dx[0]+1/dx[1]),
             D/V[1]/dx[1]*exp(c-b)],
            [0,
             D/V[2]/dx[1]*exp(b-c),
             -fref[2]/unlogt-D/dx[1]/V[2]],
        ])
    else:
        jref = np.array([
            [-D/dx[0]/V[0],             D/dx[0]/V[0],             0],
            [ D/dx[0]/V[1],  D/V[1]*(-1/dx[0]-1/dx[1]),  D/V[1]/dx[1]],
            [            0,              D/V[2]/dx[1],  -D/dx[1]/V[2]],
        ])

    if logt:
        jref *= t0

    rd.dense_jac_rmaj(t, y, jout)
    assert np.allclose(jout, jref)

    jout_bnd = np.zeros((3,N), order='F')
    rd.banded_packed_jac_cmaj(t, y, jout_bnd)
    jref_bnd = _get_banded(jref,1,N)
    assert np.allclose(jout_bnd, jref_bnd)


@pytest.mark.parametrize("geom", (FLAT, SPHERICAL, CYLINDRICAL))
def test_ReactionDiffusion__3_reactions_4_species_5_bins_k_factor(geom):
    # r[0]: A + B -> C
    # r[1]: D + C -> D + A + B
    # r[2]: B + B -> D
    pi = 3.141592653589793
    #              r[0]     r[1]    r[2]
    stoich_reac = [[0,1],   [2,3], [1,1]]
    stoich_prod = [  [2], [0,1,3],   [3]]
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

    #(r[0], r[1]) modulations over bins
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
        if bi > 0: # previous
            Cp = y0[si+n*(bi-1)]
        else:
            Cp = C # Neumann boundary conditions

        if bi < N-1:
            Cn = y0[si+n*(bi+1)]
        else:
            Cn = C # Neumann boundary conditions

        f = 0.0
        if bi > 0:
            f -= D[si]*A[bi]*(C-Cp)/dx[bi-1]
        if bi < N-1:
            f += D[si]*A[bi+1]*(Cn-C)/dx[bi]
        return f

    r = [
        [k[0]*bin_k_factor[bi][0]*C_(0,bi)*C_(1,bi) for\
         bi in range(N)],
        [k[1]*bin_k_factor[bi][1]*C_(3,bi)*C_(2,bi) for \
         bi in range(N)],
        [k[2]*C_(1,bi)**2 for bi in range(N)],
    ]

    ref_f = np.array([
        [
           -r[0][bi] + r[1][bi] +flux(0, bi)/V[bi],
           -r[0][bi] + r[1][bi] -2*r[2][bi] + flux(1, bi)/V[bi],
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
            totl = stoich_prod[ri].count(lri) - \
                   stoich_reac[ri].count(lri)
            if totl == 0: continue
            actv = stoich_actv[ri].count(lci)
            if actv == 0: continue
            v += actv*totl*r[ri][bi]/C_(lci, bi)
        return v

    def jac_elem(ri, ci):
        bri, bci = ri // n, ci // n
        lri, lci = ri  % n, ci  % n
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

    ref_banded_j = _get_banded(ref_j, n, N)

    jout_bnd_packed_cmaj = np.zeros((2*n+1, n*N), order='F')
    rd.banded_packed_jac_cmaj(0.0, y0, jout_bnd_packed_cmaj)

    plot = False
    if plot:
        import matplotlib.pyplot as plt
        from chemreac.util import coloured_spy
        fig = plt.figure()
        ax = fig.add_subplot(3,1,1)
        coloured_spy(ref_banded_j, ax=ax)
        plt.title('ref_banded_j')
        ax = fig.add_subplot(3,1,2)
        coloured_spy(jout_bnd_packed_cmaj, ax=ax)
        plt.title('jout_bnd_packed_cmaj')
        ax = fig.add_subplot(3,1,3)
        coloured_spy(ref_banded_j-jout_bnd_packed_cmaj, ax=ax)
        plt.title('diff')
        plt.show()

    assert np.allclose(jout_bnd_packed_cmaj, ref_banded_j)

    jout_bnd_padded_cmaj = np.zeros((3*n+1, n*N), order='F')
    rd.banded_padded_jac_cmaj(0.0, y0, jout_bnd_padded_cmaj)
    assert np.allclose(jout_bnd_padded_cmaj[n:,:], ref_banded_j)
