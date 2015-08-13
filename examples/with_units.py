#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Demonstration of built-in unit conversion
-----------------------------------------

::

 $ python with_units.py --help

.. exec::
   echo "::\\n\\n"
   python examples/examples/with_units.py --help | sed "s/^/   /"

"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

import numpy as np

from chemreac import ReactionDiffusion
from chemreac.integrate import Integration
from chemreac.units import second, mole, metre, molar, SI_base


def main(logy=False, logt=False, unit_registry=None):
    # A -> B
    n = 2
    k0 = 0.13 * second**-1
    if unit_registry is None:
        unit_registry = SI_base
    rd = ReactionDiffusion(n, [[0]], [[1]], k=[k0], logy=logy, logt=logt,
                           unit_registry=unit_registry)
    y0 = [3.0e3 * mole/metre**3, 1.0*molar]
    t0, tend, nt = 5.0*second, 17.0*second, 42
    tout = np.linspace(t0, tend, nt+1)

    Cref = molar*np.array([3.0*np.exp(-0.13*second**-1*(tout-t0)),
                           1.0 + 3.0*(1 - np.exp(
                               -0.13*second**-1*(tout-t0)))]).transpose()

    # scipy
    integr1 = Integration('scipy', rd, y0, tout)
    return integr1, Cref, rd


if __name__ == '__main__':
    import argh
    argh.dispatch_command(main, output_file=None)
