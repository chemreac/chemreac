# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

from ._transforms import sigmoid_algebraic_4

from scipy.special import erf

def sigmoid_algebraic(x, limit=350.0, n=2):
    # IEEE 754 double precision has Emin, Emax = (-1022, 1023)
    # 1022 * log(2) = 708.39... so halfway there (~350) seems reasonable.
    return x/(1+(x/limit)**n)**(1/n)

half_sqrt_pi = 0.8862269254527580137
def sigmoid_erf(x, limit=350.0):
    return limit*erf(half_sqrt_pi*x/limit)
