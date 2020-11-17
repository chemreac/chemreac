# -*- coding: utf-8 -*-

import sys

import pytest

from aqueous_radiolysis import integrate_rd


@pytest.mark.skipif(sys.version_info[0] > 2, reason='pytest/python3 incompatible with this test')
def test_integrate_rd():
    integr = integrate_rd()
    assert integr.info['success']
