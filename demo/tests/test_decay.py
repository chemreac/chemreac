# -*- coding: utf-8 -*-

from __future__ import (
    print_function, division, absolute_import, unicode_literals
)

from math import exp

import numpy as np
import pytest

from decay import integrate_rd


def _test_rd(forgiveness=100, **kwargs):
    yout, yref, rd, info = integrate_rd(**kwargs)
    for i in range(rd.n):
        try:
            atol = info['atol'][i]
        except:
            atol = info['atol']

        try:
            rtol = info['rtol'][i]
        except:
            rtol = info['rtol']

        if rd.logy:
            ymean = np.mean(yout[..., i])
            # exp(abserr + relerr*ytrue) == exp(abserr)*exp(relerr*ytrue)
            rtol = exp(atol)
            atol = rtol*ymean
        assert np.allclose(yout[..., i], yref[..., i], rtol*forgiveness, atol*forgiveness)


def test_linear():
    _test_rd

def test_logt():
    _test_rd(logt=True)

def test_logy_logt():
    _test_rd(logy=True, logt=True)

@pytest.mark.xfail
def test_logy():
    _test_rd(logy=True)

# Manually tweaking
def test_logy_tweak():
    _test_rd(logy=True, small=5)
