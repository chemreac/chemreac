# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import, unicode_literals

from itertools import product

import numpy as np
import pytest

from analytic_diffusion import integrate_rd

tf = (True, False)
tol = {3: 1e5, 5: 1e4, 7: 1e2} # determined from analytic_N_scaling demo
@pytest.mark.parametrize('params', list(product('fcs', tf, tf, tf, [0, .1], [3, 5, 7])))
def test_gaussian_diffusion(params):
    g, ly, lt, r, k, ns = params
    res = integrate_rd(geom=g, logt=lt, logy=ly, N=128, random=r, k=k, nstencil=ns,
                       atol=1e-6, rtol=1e-6)
    for ave_rmsd_over_atol in res[3]:
        if r:
            forgiveness = 5 if ns < 7 else 50
            # Randomized grid has lower convergence order.
            for i in range(6):
                if np.all(ave_rmsd_over_atol > tol[ns]*forgiveness):
                    ave_rmsd_over_atol = integrate_rd(
                        geom=g, logt=lt, logy=ly, N=128, random=r, k=k, nstencil=ns,
                        atol=1e-6, rtol=1e-6)[3]
                else:
                    break
            assert np.all(ave_rmsd_over_atol < tol[ns]*forgiveness)
        else:
            assert np.all(ave_rmsd_over_atol < tol[ns])
