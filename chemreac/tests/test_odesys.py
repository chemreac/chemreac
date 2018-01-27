# -*- coding: utf-8 -*-
import numpy as np
from chemreac import ReactionDiffusion
from chemreac.odesys import ODESys


def test_decay():
    # A -> B
    n = 2
    k0 = 0.13
    names = ['A', 'B']
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], substance_names=names)
    odesys = rd.as_odesys()

    y0 = [3.0, 1.0]
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt+1)

    result = odesys.integrate(tout, dict(zip(names, y0)))
    yref = np.array([y0[0]*np.exp(-k0*(tout-t0)),
                     y0[1]+y0[0]*(1-np.exp(-k0*(tout-t0)))]).transpose()
    assert np.allclose(result.yout, yref)
