#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This examples demonstrates how to use different libraries for performing
the numerical integration (stepping with error control). The following
integrators are used:

- VODE (via SciPy)
- CVODE (Sundials called directly)
- odeint (via pyodeint)
- gsl (via pygslodeiv2)


"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from chemreac import ReactionDiffusion, DENSE
from chemreac.integrate import Integration


def main(logy=False, logt=False, forgive_fact=1):
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

    print('scipy')
    print('-----')
    integr1 = Integration('scipy', rd, np.asarray(y0), np.asarray(tout))
    print(integr1.info)
    assert np.allclose(integr1.Cout[:, 0, :], Cref(tout))
    print()

    print('sundials')
    print('--------')
    integr2 = Integration('sundials', rd, np.asarray(y0), np.asarray(tout),
                          atol=[1e-8, 1e-8], rtol=1e-8, method='bdf')
    assert np.allclose(integr2.Cout[:, 0, :], Cref(tout))
    print()

    # rk4 - fixed step size with 42 steps will give poor accuracy
    print('rk4 - fixed step')
    print('----------------')
    integr3 = Integration('rk4', rd, np.asarray(y0), np.asarray(tout))
    if logt:
        assert np.allclose(integr3.Cout[:, 0, :], Cref(tout),
                           atol=4e-2, rtol=4e-2)
    else:
        assert np.allclose(integr3.Cout[:, 0, :], Cref(tout),
                           atol=5e-3, rtol=5e-3)
    print()

    print('odeint')
    print('------')
    for method, forgive in [('dopri5', 0.9), ('bs', 0.4),
                            ('rosenbrock4', 0.3)]:
        forgive *= forgive_fact
        print(method, forgive)
        atol, rtol = 1e-8, 1e-8
        integr4 = Integration('pyodeint', rd, np.asarray(y0), (t0, tend),
                              method=method, mode=DENSE, dense_output=True,
                              atol=atol, rtol=rtol)
        assert np.allclose(integr4.Cout[:, 0, :], Cref(integr4.tout),
                           atol=atol*forgive, rtol=rtol*forgive)
    print()

    print('gslodeiv2')
    print('---------')
    for method, forgive in [('bsimp', 1e-5), ('msadams', 9), ('rkf45', 0.5),
                            ('rkck', 0.3), ('rk8pd', 0.3), ('rk4imp', 1.3),
                            ('msbdf', 25)]:
        if method == 'rk4imp' and logt and not logy:
            continue  # known not to work
        forgive *= forgive_fact
        print(method, forgive)
        atol, rtol = 1e-8, 1e-8
        integr5 = Integration('pygslodeiv2', rd, np.asarray(y0), (t0, tend),
                              mode=DENSE, dense_output=True, atol=atol,
                              rtol=rtol, method=method)
        assert np.allclose(integr5.Cout[:, 0, :], Cref(integr5.tout),
                           atol=atol*forgive, rtol=rtol*forgive)


if __name__ == '__main__':
    import argh
    argh.dispatch_command(main)
