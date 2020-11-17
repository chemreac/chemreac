# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

from itertools import product

import numpy as np
import pytest

from steady_state import integrate_rd


TR_FLS = (True, False)


@pytest.mark.parametrize('params', list(product(TR_FLS, TR_FLS, TR_FLS,
                                                # TR_FLS, TR_FLS,
                                                [3], [1e-4])))
def test_steady_state(params):
    # ly, lt, r, lr, rr, ns, forgiveness = params
    ly, lt, r, ns, forgiveness = params
    lr, rr = False, False
    # notice forgiveness << 1
    # this is because static conditions is very simple to integrate
    atol = 1e-6
    tout, yout, info, ave_rmsd_over_atol, rd = integrate_rd(
        geom='f', logt=lt, logy=ly, N=128, random=r, nstencil=ns,
        lrefl=lr, rrefl=rr, atol=atol, rtol=1e-6)
    assert info['success']
    for rmsd in ave_rmsd_over_atol:
        assert np.all(rmsd < forgiveness)


@pytest.mark.parametrize('params', list(product(
    TR_FLS, TR_FLS, TR_FLS, [5, 7])))
def test_steady_state__high_stencil(params):
    ly, lt, r, nstencil = params
    test_steady_state((ly, lt, r,  # False, False,
                       nstencil, 3e-2))


# @pytest.mark.xfail
# @pytest.mark.parametrize('params', list(product(TR_FLS, TR_FLS,
#     TR_FLS, 'cs')))
# def test_steady_state__curved_geom(params):
#     ly, lt, r, geom = params
#     atol = 1e-6
#     tout, yout, info, ave_rmsd_over_atol, rd = integrate_rd(
#         geom=geom, logt=lt, logy=ly, N=128, random=r,
#         lrefl=True, rrefl=True, atol=atol, rtol=1e-6)
#     assert info['success']
#     for rmsd in ave_rmsd_over_atol:
#         assert np.all(rmsd < 1.0)
