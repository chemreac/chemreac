# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

from itertools import product

import numpy as np
import pytest

from analytic_diffusion import integrate_rd
from chemreac.util.testing import veryslow

TR_FLS = (True, False)

COMBOS = list(product('fcs', TR_FLS, TR_FLS, [False], [.1], [3]))
EXTRA_COMBOS = list(product('fcs', TR_FLS, TR_FLS, TR_FLS, [0, .2], [5, 7]))

tol = {3: 1e5, 5: 1e4, 7: 1e2}  # determined from analytic_N_scaling demo


def _test_gaussian_diffusion(params):
    g, ly, lt, r, k, ns = params
    tout, yout, info, ave_rmsd_over_atol, sys = integrate_rd(
        geom=g, logt=lt, logy=ly, N=128, random=r, k=k, nstencil=ns,
        atol=1e-6, rtol=1e-6)
    assert info['success']
    for rmsd in ave_rmsd_over_atol:
        if r:
            forgiveness = 5 if ns < 7 else 50
            # Randomized grid has lower convergence order.
            for i in range(6):
                if np.all(rmsd > tol[ns]*forgiveness):
                    rmsd = integrate_rd(
                        geom=g, logt=lt, logy=ly, N=128, random=r, k=k,
                        nstencil=ns, atol=1e-6, rtol=1e-6)[3]
                else:
                    break
            assert np.all(rmsd < tol[ns]*forgiveness)
        else:
            assert np.all(rmsd < tol[ns])


@pytest.mark.parametrize('params', COMBOS)
def test_gaussian_diffusion(params):
    _test_gaussian_diffusion(params)


@veryslow
@pytest.mark.parametrize('params', EXTRA_COMBOS)
def test_gaussian_diffusion_extended(params):
    _test_gaussian_diffusion(params)
