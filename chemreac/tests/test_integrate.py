#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

try:
    import cPickle as pickle
except ImportError:
    import pickle

from itertools import product
import os

import numpy as np
import pytest

from chemreac import ReactionDiffusion
from chemreac.integrate import run, Integration
from chemreac.util.testing import veryslow
from chemreac.units import metre, molar, umol, hour, day, SI_base_registry

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

    integr = run(rd, y0, tout)
    yref = np.array([y0[0]*np.exp(-k0*(tout-t0)),
                     y0[1]+y0[0]*(1-np.exp(-k0*(tout-t0)))]).transpose()
    assert np.allclose(integr.Cout[:, 0, :], yref)


def test_decay_solver_kwargs_env():
    key = 'CHEMREAC_INTEGRATION_KWARGS'
    try:
        ori_val = os.environ.pop(key)
    except KeyError:
        ori_val = None
    os.environ[key] = "{'rtol': 1e-9}"
    test_decay(LOG_COMOBS[0])
    if ori_val is not None:
        os.environ[key] = ori_val
    else:
        os.environ.pop(key)


def test_autodimerization():
    # A + A -> B
    from chemreac.chemistry import (
        Reaction, ReactionSystem, mk_sn_dict_from_names
    )
    sbstncs = mk_sn_dict_from_names('AB')
    k = 3.0
    r1 = Reaction({'A': 2}, {'B': 1}, k)
    rsys = ReactionSystem([r1], sbstncs)
    rd = ReactionDiffusion.from_ReactionSystem(rsys)
    t = np.linspace(0, 5, 3)
    A0, B0 = 1.0, 0.0
    integr = run(rd, [A0, B0], t)
    Aref = 1/(1/A0+2*k*t)
    yref = np.vstack((Aref, (A0-Aref)/2)).transpose()
    assert np.allclose(integr.yout[:, 0, :], yref)


@veryslow
@pytest.mark.parametrize("log_geom", product(
    LOG_COMOBS, 'fcs'))
def test_ReactionDiffusion_fields_and_g_values(log_geom):
    # modulation in x means x_center
    # A -> B # mod1 (x**2)
    # C -> D # mod1 (x**2)
    # E -> F # mod2 (sqrt(x))
    # G -> H # no modulation
    (logy, logt), geom = log_geom
    k = np.array([3.0, 7.0, 13.0, 22.0])
    N = 5
    n = 8
    D = np.zeros(n)
    x = np.linspace(3, 7, N+1)
    xc = x[:-1] + np.diff(x)/2
    fields = [[xc[i]*xc[i] for i in range(N)],
              [xc[i]*xc[i] for i in range(N)],
              [xc[i]**0.5 for i in range(N)]]
    g_values = [[-k[0], k[0], 0, 0, 0, 0, 0, 0],
                [0, 0, -k[1], k[1], 0, 0, 0, 0],
                [0, 0, 0, 0, -k[2], k[2], 0, 0]]
    g_value_parents = [0, 2, 4]
    rd = ReactionDiffusion(
        n,
        [[i] for i in range(n-2, n, 2)],
        [[i] for i in range(n-1, n, 2)],
        k=k[3:], N=N, D=D, x=x, fields=fields,
        g_values=g_values, g_value_parents=g_value_parents,
        logy=logy, logt=logt, geom=geom
    )
    assert rd.n == n
    assert rd.N == N
    assert np.allclose(rd.D, D)
    assert np.allclose(rd.k, k[3:])
    assert np.allclose(rd.xcenters, xc)
    assert np.allclose(rd.fields, fields)
    assert np.allclose(rd.g_values, g_values)
    y0 = np.array([[13.0, 23.0, 32.0, 43.0, 12.0, 9.5, 17.0, 27.5]*N])
    y0 = y0.flatten()
    t0, tend, nt = 1.0, 1.1, 42
    tout = np.linspace(t0, tend, nt+1)

    integr = run(rd, y0, tout, atol=1e-11, rtol=1e-11, with_jacobian=True)

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
    assert np.allclose(integr.Cout.flatten(), yref.flatten())


@pytest.mark.parametrize("N_wjac_geom", product(
    [64, 128], [False, True], 'fcs'))
def test_integrate__only_1_species_diffusion__mass_conservation(N_wjac_geom):
    N, wjac, geom = N_wjac_geom
    # Test that mass convervation is fulfilled wrt diffusion.
    x = np.linspace(0.01*N, N, N+1)
    y0 = (x[0]/2/N+x[1:]/N)**2

    sys = ReactionDiffusion(1, [], [], [], N=N, D=[0.02*N], x=x, geom=geom,
                            nstencil=3, lrefl=True, rrefl=True)

    tout = np.linspace(0, 10.0, 50)
    atol, rtol = 1e-6, 1e-8
    integr = run(sys, y0, tout, atol=atol, rtol=rtol, with_jacobian=wjac,
                 method='adams')
    yout = integr.yout[:, :, 0]
    x /= N
    if geom == 'f':
        yprim = yout*(x[1:]**1 - x[:-1]**1)
    elif geom == 'c':
        yprim = yout*(x[1:]**2 - x[:-1]**2)
    elif geom == 's':
        yprim = yout*(x[1:]**3 - x[:-1]**3)
    else:
        raise

    ybis = np.sum(yprim, axis=1)

    assert np.allclose(np.average(ybis), ybis, atol=atol, rtol=rtol)


@pytest.mark.parametrize("log", LOG_COMOBS)
def test_integrators(log):
    logy, logt = log
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt+1)

    # Update the dict if more integrators are added:
    solver_kwargs = {
        'scipy1': {
            'atol': [1e-8, 1e-8],
            'rtol': 1e-8,
            'tout': tout
        },
        'scipy2': {
            'atol': 1e-8,
            'rtol': 1e-8,
            'tout': (t0, tend),
            'dense_output': True
        },
        'cvode1': {
            'atol': [1e-8, 1e-8],
            'rtol': 1e-8,
            'method': 'bdf',
            'tout': tout
        },
        'cvode2': {
            'atol': 1e-8,
            'rtol': 1e-8,
            'method': 'adams',
            'tout': tout,
            'C0_is_log': True
        }
    }

    # A -> B
    n = 2
    k0 = 0.13
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], logy=logy, logt=logt)
    y0 = [3.0, 1.0]

    results = []
    for solver, kwargs in solver_kwargs.items():
        _y0 = np.log(y0) if kwargs.get('C0_is_log', False) else y0
        integr = Integration(rd, _y0, integrator=solver[:-1], **kwargs)
        if not kwargs.get('dense_output', False):
            results.append(integr.Cout)

    for result in results[1:]:
        assert np.allclose(results[0][0], result[0])


def test_pickle_Integration():
    # A -> B
    n = 2
    k0 = 0.13
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0])
    y0 = [3.0, 1.0]
    t0, tend, nt = 5.0, 17.0, 42
    integr = Integration(rd, y0, tout=np.linspace(t0, tend, nt+1),
                         integrator='scipy')
    Cout = integr.Cout.copy()
    s = pickle.dumps(integr)
    integr2 = pickle.loads(s)
    assert np.allclose(integr2.Cout, Cout)


def test_integrate_nondimensionalisation():
    # 2A -> B
    rd = ReactionDiffusion.nondimensionalisation(
        2, [[0, 0]], [[1]], [2e-9/(umol/metre**3)/hour],
        unit_registry=SI_base_registry)
    C0 = [3*molar, 4*molar]
    tout = np.linspace(0, 1)*day
    integr = Integration.nondimensionalisation(
        rd, C0, tout, integrator='scipy')

    k_m3_p_mol_p_sec = 2e-3/3600
    t_sec = np.linspace(0, 24*3600)
    C0_mol_p_m3 = [3000, 4000]
    Cref_mol_p_m3 = np.empty(integr.Cout.squeeze().shape)
    Cref_mol_p_m3[:, 0] = 1/(C0_mol_p_m3[0]**-1 + 2*k_m3_p_mol_p_sec*t_sec)
    missing_A = (C0_mol_p_m3[0] - Cref_mol_p_m3[:, 0])
    Cref_mol_p_m3[:, 1] = C0_mol_p_m3[1] + missing_A/2
    assert np.allclose(integr.tout, t_sec)
    print(integr.Cout.squeeze())
    print(Cref_mol_p_m3)
    print(integr.Cout.squeeze() - Cref_mol_p_m3)
    assert np.allclose(integr.Cout.squeeze(), Cref_mol_p_m3)
