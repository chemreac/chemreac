#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np
import pytest

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL


@pytest.mark.xfail
def test_ReactionDiffusion__f__wrong_fout_dimension():
    y0 = np.array([2.0, 3.0])
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k])
    fout = np.ones((1,))*99 # fout too small
    rd.f(0.0, y0, fout)


@pytest.mark.parametrize("N", range(1,5))
def test_ReactionDiffusion__only_1_reaction(N):
    y0 = np.array([2.0, 3.0]*N)
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, [[0]], [[1]], [k], N, D=[0.0, 0.0])
    fout = np.ones((2*N,))*99
    rd.f(0.0, y0, fout)

    for i in range(N):
        assert np.allclose(fout[i*2:(i+1)*2], np.array([-10.0, 10.0]))

@pytest.mark.parametrize("N", range(1,5))
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


def test_ReactionDiffusion__only_1_species_diffusion():
    # Diffusion without reactionsn
    D = 17.0
    y0 = np.array([23.0, 27.0])
    x = [5.0, 7.0, 13.0]
    rd = ReactionDiffusion(1, [], [], [], D=[D], x=x)
    fout = np.ones((2,))*99
    rd.f(0.0, y0, fout)
    J = D*(y0[0]-y0[1])/(0.5*(x[2]-x[0]))
    fref = np.array([-J/(x[1]-x[0]), J/(x[2]-x[1])])
    assert np.allclose(fout, fref)

def test_ReactionDiffusion__actv():
    pass

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
    y0 = np.array([2.5, 1.2, 3.2, 4.3,
                   2.7, 0.8, 1.6, 2.4,
                   3.1, 0.3, 1.5, 1.8,
                   3.3, 0.6, 1.6, 1.4,
                   3.6, 0.9, 1.7, 1.2])
    x = np.array([11.0, 13.0, 17.0, 23.0, 29.0, 37.0])
    k = [31.0, 37.0, 41.0]

    #(r[0], r[1]) modulations
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

    jout_bnd_cmaj = np.zeros((n*N, n*N), order='F')
    jout_pck_bnd_cmaj = np.zeros((n*N, n*N), order='F')
