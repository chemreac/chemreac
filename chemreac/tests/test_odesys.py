# -*- coding: utf-8 -*-
import numpy as np
from chemreac import ReactionDiffusion

analytic = [
    lambda y0, k, t: (
        y0[0] * np.exp(-k[0]*t)),
    lambda y0, k, t: (
        y0[1] * np.exp(-k[1] * t) + y0[0] * k[0] / (k[1] - k[0]) *
        (np.exp(-k[0]*t) - np.exp(-k[1]*t))),
    lambda y0, k, t: (
        y0[2] + y0[1] * k[1] / (-k[1]) *
        (np.exp(-k[1]*t) - 1) +
        k[1] * k[0] * y0[0] / (k[1] - k[0]) *
        (1 / (-k[0]) * (np.exp(-k[0]*t) - 1) -
         1 / (-k[1]) * (np.exp(-k[1]*t) - 1)))
]


def _get_odesys(k):
    names = ['A', 'B', 'C']
    pns = ['kA', 'kB']
    rd = ReactionDiffusion(len(names), [[0], [1]], [[1], [2]], k=k,
                           substance_names=names, param_names=pns)
    return rd._as_odesys(k_from_params=lambda p: [p[k] for k in pns])


def test_decay():
    kA = 0.13
    odesys = _get_odesys([kA, 0])
    y0 = dict(A=3, B=1, C=0)
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt+1)
    result = odesys.integrate(tout, y0)
    yref = np.array([y0['A']*np.exp(-kA*(tout-t0)),
                     y0['B']+y0['A']*(1-np.exp(-kA*(tout-t0)))]).transpose()
    assert np.allclose(result.yout[:, :2], yref)


def test_decay_params():
    odesys = _get_odesys([0, 0])
    y0 = 42, 7, 4
    k = .7, .3
    ic = dict(zip(odesys.names, y0))
    p = dict(zip('kA kB'.split(), k))
    tout, yout, info = odesys.integrate([0, 5], ic, p)
    yref = np.array([a(y0, k, tout) for a in analytic]).transpose()
    assert np.allclose(yout, yref)
