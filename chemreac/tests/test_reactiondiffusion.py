#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from chemreac import PyReactionDiffusion as ReactionDiffusion


def test_ReactionDiffusion__only_1_reaction():
    y0 = np.array([2.0, 3.0])
    k = 5.0
    # A -> B
    rd = ReactionDiffusion(2, 1, [[0]], [[1]], [k])
    fout = np.ones((2,))*99
    rd.f(0.0, y0, fout)
    assert np.allclose(fout, np.array([-10.0, 10.0]))


def test_ReactionDiffusion__only_1_species_diffusion():
    # Diffusion without reactionsn
    D = 17.0
    y0 = np.array([23.0, 27.0])
    x = [5.0, 7.0, 13.0]
    rd = ReactionDiffusion(1, 2, [], [], [], [D], x=x)
    fout = np.ones((2,))*99
    rd.f(0.0, y0, fout)
    print(fout)
    J = D*(y0[0]-y0[1])/(0.5*(x[2]-x[0]))
    fref = np.array([-J/(x[1]-x[0]), J/(x[2]-x[1])])
    assert np.allclose(fout, fref)


if __name__ == '__main__':
    test_ReactionDiffusion__only_1_reaction()
    test_ReactionDiffusion__only_1_species_diffusion()
