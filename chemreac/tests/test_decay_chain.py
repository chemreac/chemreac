# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import

import numpy as np

from chemreac import ReactionDiffusion
from .test_reactiondiffusion import _test_f, _test_dense_jac_rmaj


def test_dc():
    n = 7
    k = [1]*(n-1)
    rd = ReactionDiffusion(n, [[i] for i in range(n-1)],
                           [[i] for i in range(1, n)], k=k)
    fout = rd.alloc_fout()
    y0 = [0]*n
    y0[0] = 1
    y_ = [0] + y0
    y0 = np.asarray(y0, dtype=np.float64)
    rd.f(0, y0, fout)
    pos, neg = [0]+k, k+[0]
    _test_f(rd, 0, y0)
    fref = [pos[i]*y_[i] - neg[i]*y_[i+1] for i in range(n)]
    assert np.allclose(fref, fout)

    jout = rd.alloc_jout(order='C')
    rd.dense_jac_rmaj(0, y0, jout)
    jref = np.zeros_like(jout)
    for i in range(n):
        if i < n - 1:
            jref[i, i] = -k[i]
        if i > 0:
            jref[i, i-1] = k[i-1]
    assert np.allclose(jref, jout)
    _test_dense_jac_rmaj(rd, 0, y0)
