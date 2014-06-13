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
    kwargs = {'geom': geom, 'logt': lt, 'logy': ly, 'logx': lx,
              'random': r, 'nstencil': ns, 'lrefl': lr, 'rrefl': rr}

    # zero offset
    tout, yout, info, sys = integrate_rd(D=1e-3, N=64,
                                         offset=0, **kwargs)
    if ly:
        delta_y = np.exp(yout[:, :, 0]) - np.exp(yout[:, :, 1])
    else:
        delta_y = yout[:, :, 0] - yout[:, :, 1]
    assert np.all(np.abs(delta_y) < 1e-3)

    # separated gaussians
    tout, yout, info, sys = integrate_rd(
        D=0.0, t0=1e-6, tend=7, x0=0.1, xend=1.0, N=64,
        offset=.25, mobility=3e-7, nt=25, sigma_q=101, **kwargs)
    mass_consv = np.sum(np.exp(yout) if ly else yout, axis=1)
    #assert np.all(mass_consv - mass_consv[[0],:] < 1e-6)
    if geom == 'f':
        efield_i = sys.calc_efield(yout[0,:,:].flatten())
        efield_f = sys.calc_efield(yout[-1,:,:].flatten())
        assert np.max(np.abs(efield_i)) - np.max(np.abs(efield_f)) < 1e-6
