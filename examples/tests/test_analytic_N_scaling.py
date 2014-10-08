# -*- coding: utf-8 -*-

from analytic_N_scaling import main
from chemreac.util.testing import veryslow


@veryslow
def test_integrate_rd():
    main(nNs=5)
