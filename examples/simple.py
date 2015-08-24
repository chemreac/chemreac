#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

import numpy as np

from chemreac import ReactionDiffusion, DENSE
from chemreac.integrate import Integration


def main(logy=False, logt=False):
    # A -> B
    n = 2
    k0 = 0.13
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], logy=logy, logt=logt)
    y0 = [3.0, 1.0]
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt)

    def Cref(tarr):
        return np.array([
            y0[0]*np.exp(-k0*(tarr-tarr[0])),
            y0[1] + y0[0]*(1 - np.exp(-k0*(tarr-tarr[0])))]).transpose()

    # scipy
    integr1 = Integration('scipy', rd, np.asarray(y0), np.asarray(tout))
    assert np.allclose(integr1.Cout[:, 0, :], Cref(tout))

    # sundials
    integr2 = Integration('sundials', rd, np.asarray(y0), np.asarray(tout),
                          atol=[1e-8, 1e-8], rtol=1e-8, method='bdf')
    assert np.allclose(integr2.Cout[:, 0, :], Cref(tout))

    # rk4 - fixed step size with 42 steps will give poor accuracy
    integr3 = Integration('rk4', rd, np.asarray(y0), np.asarray(tout))
    if logt:
        assert np.allclose(integr3.Cout[:, 0, :], Cref(tout),
                           atol=4e-2, rtol=4e-2)
    else:
        assert np.allclose(integr3.Cout[:, 0, :], Cref(tout),
                           atol=5e-3, rtol=5e-3)

    # odeint
    integr4 = Integration('odeint', rd, np.asarray(y0), (t0, tend),
                          mode=DENSE, dense_output=True, atol=1e-9, rtol=1e-9)
    assert np.allclose(integr4.Cout[:, 0, :], Cref(integr4.tout),
                       atol=1e-5, rtol=1e-5)


if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
