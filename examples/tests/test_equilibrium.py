# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *


from math import exp

import pytest

import equilibrium
from chemreac.util.testing import _test_rd_integration_run


def test_analytic_x():
    assert abs(equilibrium.analytic_x(3, 5, 13, 7, 2, 11) -
               1.18768491286119571) < 1e-14


def test_linear():
    _test_rd_integration_run(equilibrium.integrate_rd)


def test_logt():
    _test_rd_integration_run(equilibrium.integrate_rd, 10, logt=True)


def test_logy_logt():
    _test_rd_integration_run(equilibrium.integrate_rd, 40, rtol='1e-10',
                             logy=True, logt=True)
    _test_rd_integration_run(equilibrium.integrate_rd, 20, tend=2.0,
                             atol='1e-6', rtol='1e-13', logy=True, logt=True)


def test_logy():
    _test_rd_integration_run(equilibrium.integrate_rd, 30, logy=True)
