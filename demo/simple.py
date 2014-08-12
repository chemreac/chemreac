#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

import argh
import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import run


def main(logy=False, logt=False):
    # A -> B
    n = 2
    k0 = 0.13
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], logy=logy, logt=logt)
    y0 = [3.0, 1.0]
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt+1)

    # y = np.log(y0) if logy else y0
    # t = np.log(tout) if logt else tout
    # yout, info = run(rd, y, t)
    # yout = np.exp(yout) if logy else yout

    # yref = np.array([y0[0]*np.exp(-k0*(tout-t0)),
    #                  y0[1]+y0[0]*(1-np.exp(-k0*(tout-t0)))]).transpose()
    # assert np.allclose(yout[:, 0, :], yref)

    # sundials
    from chemreac._chemreac import sundials_direct
    print(rd.auto_efield)
    yout2 = sundials_direct(rd, [1e-8, 1e-8], 1e-8, 'bdf', np.asarray(y0), tout)
    assert np.allclose(yout2[:, 0, :], yref)


if __name__ == '__main__':
    argh.dispatch_command(main)
