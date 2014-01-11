#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pytest

from chemreac import PyReactionDiffusion as ReactionDiffusion

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
