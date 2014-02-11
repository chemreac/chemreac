#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import, unicode_literals


import os
from itertools import product

import numpy as np
import pytest

from chemreac import ReactionDiffusion, FLAT, SPHERICAL, CYLINDRICAL
from chemreac.integrate import run
from chemreac.serialization import load

"""
Tests the integration routine for the
chemical reaction system. (no diffusion)
"""

LOG_COMOBS = list(product([True, False], [True, False]))

@pytest.mark.parametrize("log", LOG_COMOBS)
def test_decay(log):
    # A -> B
    n = 2
    logy, logt = log
    k0 = 0.13
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], logy=logy, logt=logt)
    y0 = [3.0, 1.0]
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt+1)

    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(rd, y, t)
    if logy: yout = np.exp(yout)

    yref = np.array([y0[0]*np.exp(-k0*(tout-t0)),
                     y0[1]+y0[0]*(1-np.exp(-k0*(tout-t0)))]).transpose()
    assert np.allclose(yout, yref)


@pytest.mark.parametrize("log_geom", product(LOG_COMOBS, (FLAT, SPHERICAL, CYLINDRICAL)))
def test_ReactionDiffusion__bin_k_factor(log_geom):
    # A -> B # mod1 (x**2)
    # C -> D # mod1 (x**2)
    # E -> F # mod2 (sqrt(x))
    # G -> H # no modulation
    (logy, logt), geom = log_geom
    k = np.array([3.0, 7.0, 13.0, 22.0])
    N = 5
    n = 8
    nr = 4
    D = np.zeros(n)
    x = np.linspace(3,7,N)
    bkf = [(x[i]*x[i], x[i]**0.5) for i in range(N)]
    bkf_span = [2,1]
    rd = ReactionDiffusion(
        n,
        [[i] for i in range(0,n,2)],
        [[i] for i in range(1,n,2)],
        k=k, N=N, D=D, bin_k_factor=bkf,
        bin_k_factor_span=bkf_span,
        logy=logy, logt=logt, geom=geom
    )
    y0 = np.array([[13.0, 23.0, 32.0, 43.0, 12.0, 9.5, 17.0, 27.5]*N]).flatten()
    t0, tend, nt = 1.0, 1.1, 42
    tout = np.linspace(t0, tend, nt+1)

    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(rd, y, t)
    if logy: yout = np.exp(yout)

    def _get_bkf(bi, ri):
        if ri < 2:
            return x[bi]*x[bi]
        elif ri < 3:
            return x[bi]**0.5
        else:
            return 1.0

    yref = np.hstack([
        np.hstack([
            np.array([
                y0[i]*np.exp(-_get_bkf(bi, i/2)*k[i/2]*(tout-t0)),
                y0[i+1]+y0[i]*(1-np.exp(-_get_bkf(bi, i/2)*k[i/2]*(tout-t0)))
            ]).transpose() for i in range(0,n,2)
        ]) for bi in range(N)])
    assert np.allclose(yout, yref)


@pytest.mark.parametrize("N", range(2,17))
def test_integrate__only_1_species_diffusion__mass_conservation(N):
    # Test that mass convervation is fulfilled wrt diffusion.
    x = np.linspace(0.1, 1.0, N+1)
    y0 = (x[0]/2+x[1:])**2

    geoms = (FLAT, SPHERICAL, CYLINDRICAL)

    for i, G in enumerate(geoms):
        sys = ReactionDiffusion(1, [], [], [], N=N, D=[0.02], x=x, geom=G)
        tout = np.linspace(0, 10.0, 50)
        yout, info = run(sys, y0, tout)
        if i == 0:
            yprim = yout
        elif i == 1:
            yprim = yout*(x[1:]**3-x[:-1]**3)
        else:
            yprim = yout*(x[1:]**2-x[:-1]**2)

        ybis = np.sum(yprim, axis=1)

        assert np.allclose(np.average(ybis), ybis)
