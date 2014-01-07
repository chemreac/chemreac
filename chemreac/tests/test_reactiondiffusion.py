#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from chemreac import PyReactionDiffusion as ReactionDiffusion


def test_ReactionDiffusion():

    # A -> B
    print('test1')
    rd1 = ReactionDiffusion(2, 1, [[0]], [[1]], [1.0], [1.0])
    y0 = np.array([1.0, 0.0])
    fout = np.ones((2,))*99
    rd1.f(0.0, y0, fout)
    assert np.allclose(fout, np.array([-1.0, 1.0]))


    # Diffusion without reactions
    print('test2')
    rd2 = ReactionDiffusion(1, 2, [], [], [], [1.0])
    rd2.f(0.0, y0, fout)
    print(fout)
    assert np.allclose(fout, np.array([-1.0, 1.0]))


if __name__ == '__main__':
    test_ReactionDiffusion()
