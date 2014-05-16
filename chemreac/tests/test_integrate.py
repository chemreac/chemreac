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
    yout = np.exp(yout) if logy else yout

    yref = np.array([y0[0]*np.exp(-k0*(tout-t0)),
                     y0[1]+y0[0]*(1-np.exp(-k0*(tout-t0)))]).transpose()
    assert np.allclose(yout[:, 0, :], yref)


def test_autodimerization():
    # A + A -> B
    from chemreac.chemistry import Reaction, ReactionSystem, mk_sn_dict_from_names
    sbstncs = mk_sn_dict_from_names('AB')
    k = 3.0
    r1 = Reaction({'A': 2}, {'B': 1}, k=k)
    rsys = ReactionSystem([r1])
    rd = rsys.to_ReactionDiffusion(sbstncs)
    t = np.linspace(0, 5, 3)
    A0, B0 = 1.0, 0.0
    yout, info = run(rd, [A0, B0], t)
    Aref = 1/(1/A0+2*k*t)
    yref = np.vstack((Aref, (A0-Aref)/2)).transpose()
    assert np.allclose(yout[:, 0, :], yref)


@pytest.mark.parametrize("log_geom", product(LOG_COMOBS, (FLAT, SPHERICAL, CYLINDRICAL)))
def test_ReactionDiffusion__bin_k_factor(log_geom):
    # modulation in x means x_center
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
    x = np.linspace(3, 7, N+1)
    xc = x[:-1] + np.diff(x)/2
    bkf = [(xc[i]*xc[i], xc[i]**0.5) for i in range(N)]
    bkf_span = [2,1]
    rd = ReactionDiffusion(
        n,
        [[i] for i in range(0,n,2)],
        [[i] for i in range(1,n,2)],
        k=k, N=N, D=D, x=x, bin_k_factor=bkf,
        bin_k_factor_span=bkf_span,
        logy=logy, logt=logt, geom=geom
    )
    assert np.allclose(rd.xcenters, xc)
    y0 = np.array([[13.0, 23.0, 32.0, 43.0, 12.0, 9.5, 17.0, 27.5]*N]).flatten()
    t0, tend, nt = 1.0, 1.1, 42
    tout = np.linspace(t0, tend, nt+1)

    y = np.log(y0) if logy else y0
    t = np.log(tout) if logt else tout
    yout, info = run(rd, y, t, atol=1e-11, rtol=1e-11)
    if logy: yout = np.exp(yout)

    def _get_bkf(bi, ri):
        if ri < 2:
            return xc[bi]*xc[bi]
        elif ri < 3:
            return xc[bi]**0.5
        else:
            return 1.0

    yref = np.hstack([
        np.hstack([
            np.array([
                y0[i]*np.exp(-_get_bkf(bi, i/2)*k[i/2]*(tout-t0)),
                y0[i+1]+y0[i]*(1-np.exp(-_get_bkf(bi, i/2)*k[i/2]*(tout-t0)))
            ]).transpose() for i in range(0, n, 2)
        ]) for bi in range(N)])
    assert np.allclose(yout.flatten(), yref.flatten())


@pytest.mark.parametrize("N_wjac_geom", product(
    [64, 128], [False, True], (FLAT, CYLINDRICAL, SPHERICAL)))
def test_integrate__only_1_species_diffusion__mass_conservation(N_wjac_geom):
    N, wjac, geom = N_wjac_geom
    # Test that mass convervation is fulfilled wrt diffusion.
    x = np.linspace(0.01*N, N, N+1)
    y0 = (x[0]/2/N+x[1:]/N)**2

    sys = ReactionDiffusion(1, [], [], [], N=N, D=[0.02*N], x=x, geom=geom, nstencil=3,
                            lrefl=True, rrefl=True)
    debug = False
    if debug: # From debugging / test design
        import matplotlib.pyplot as plt
        plt.subplot(2,1,1)
        plt.plot(sys.xcenters, y0)
        plt.xlabel('x')
        plt.ylabel('y0')

    tout = np.linspace(0, 10.0, 50)
    atol, rtol = 1e-6, 1e-8
    yout, info = run(sys, y0, tout, atol=atol, rtol=rtol, with_jacobian=wjac, method='adams')
    yout = yout[:,:,0]
    x /= N
    if geom == FLAT:
        yprim = yout*(x[1:]**1 - x[:-1]**1)
    elif geom == CYLINDRICAL:
        yprim = yout*(x[1:]**2 - x[:-1]**2)
    elif geom == SPHERICAL:
        yprim = yout*(x[1:]**3 - x[:-1]**3)
    else:
        raise

    ybis = np.sum(yprim, axis=1)

    if debug: # From debugging / test design
        plt.subplot(2,1,2)
        plt.plot(ybis-ybis[0])
        plt.xlabel('t')
        plt.ylabel('tot y - tot y0')
        plt.show()
        print(y0)
        print(np.average(ybis) - ybis)
        print(np.average(ybis))
    assert np.allclose(np.average(ybis), ybis, atol=atol, rtol=rtol)
