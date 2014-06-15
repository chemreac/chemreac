# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import numpy as np

from chemreac.util.transforms import (
    sigmoid_algebraic_4, sigmoid_algebraic, sigmoid_erf
)


def test_sigmoid_algebraic():
    assert abs(5 - sigmoid_algebraic(5, n=2)) < 1e-3
    assert abs(5 - sigmoid_algebraic(5, n=4)) < 1e-7
    assert abs(9 - sigmoid_algebraic(9, n=6)) < 1e-9
    assert abs(5 + sigmoid_algebraic(-5, n=2)) < 1e-3
    assert abs(5 + sigmoid_algebraic(-5, n=4)) < 1e-7
    assert abs(9 + sigmoid_algebraic(-9, n=6)) < 1e-9
    assert np.allclose(sigmoid_algebraic(np.array(
        [-9, -5, 5, 9]), n=6), [-9, -5, 5, 9], atol=1e-7)


def test_sigmoid_erf():
    assert abs(5 - sigmoid_erf(5)) < 1e-3
    assert abs(5 + sigmoid_erf(-5)) < 1e-3
    assert np.allclose(sigmoid_erf(np.array(
        [-5, 5])), [-5, 5], atol=1e-3)


def test_sigmoid_algebraic_4():
    assert abs(5 - sigmoid_algebraic_4(5)) < 1e-7
    assert abs(9 - sigmoid_algebraic_4(9)) < 1e-6
    assert abs(29 - sigmoid_algebraic_4(29)) < 1e-3
    assert abs(5 + sigmoid_algebraic_4(-5)) < 1e-7
    assert abs(9 + sigmoid_algebraic_4(-9)) < 1e-6
    assert abs(29 + sigmoid_algebraic_4(-29)) < 1e-3
    assert np.allclose(sigmoid_algebraic_4(np.array(
        [-9, -5, 5, 9]).reshape((2,2))), [[-9, -5], [5, 9]], atol=1e-6)
