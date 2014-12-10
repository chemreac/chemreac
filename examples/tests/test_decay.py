# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *


from math import exp

import pytest

import decay
from chemreac.util.testing import _test_rd_integration_run


def test_linear():
    _test_rd_integration_run(decay.integrate_rd)


def test_logt():
    _test_rd_integration_run(decay.integrate_rd, 30, logt=True)


def test_logy_logt():
    _test_rd_integration_run(decay.integrate_rd, 500, rtol='1e-10', logy=True,
                             logt=True, small=1e-20, t0=1e-70)
    _test_rd_integration_run(decay.integrate_rd, 500, tend=4.0, rates='1.0',
                             atol='1e-6', rtol='1e-13', logy=True, logt=True,
                             small=1e-20, t0=1e-40)


# Manually tweaking
def test_logy_tweak():
    _test_rd_integration_run(decay.integrate_rd, 10, rtol='1e-11', logy=True,
                             small=1e-6)


@pytest.mark.xfail
def test_logy():
    _test_rd_integration_run(decay.integrate_rd, 300, logy=True, small=1e-20)
