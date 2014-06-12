# -*- coding: utf-8 -*-

from __future__ import (
    print_function, division, absolute_import, unicode_literals
)

from itertools import product

import numpy as np
import pytest

from auto_efield import integrate_rd

TR_FLS = (True, False)


@pytest.mark.parametrize('params', list(product(
    TR_FLS, TR_FLS, TR_FLS, TR_FLS, TR_FLS, TR_FLS, [3, 5], 'fcs')))
def test_perfect_overlap(params):
    ly, lt, lx, r, lr, rr, ns, geom = params
    tout, yout, info, sys = integrate_rd(
        geom=geom, logt=lt, logy=ly, logx=lx,
        N=64, random=r, nstencil=ns,
        lrefl=lr, rrefl=rr, atol=1e-6, rtol=1e-6, offset=0)
    if ly:
        delta_y = np.exp(yout[:, :, 0]) - np.exp(yout[:, :, 1])
    else:
        delta_y = yout[:, :, 0] - yout[:, :, 1]
    assert np.all(np.abs(delta_y) < 1e-3)
