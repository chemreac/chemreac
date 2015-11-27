# -*- coding: utf-8 -*-
from __future__ import print_function, division, absolute_import


import numpy as np

from chemreac.symbolic import SymRD


def test_SymRD():
    rd = SymRD(2, [[0]], [[1]], k=[5.0])

    fout = rd.alloc_fout()
    rd.f(0, [7, 3], fout)
    assert np.allclose(fout, [-5*7, 5*7])

    jout = rd.alloc_jout(banded=False)
    rd.dense_jac_rmaj(0, [7, 3], jout)
    assert np.allclose(jout, [[-5.0, 0.0], [5.0, 0.0]])
