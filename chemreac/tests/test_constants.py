# -*- coding: utf-8 -*-

from chemreac.constants import get_unitless_constant


def test_get_unitless_constant():
    vp = get_unitless_constant(None, 'vacuum_permittivity')
    assert abs(vp - 8.854187817e-12) < 1e-19
