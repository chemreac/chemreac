# -*- coding: utf-8 -*-

from collections import OrderedDict

import numpy as np

from chemreac.util.stoich import (
    get_coeff_mtx, decompose_yield_into_rate_coeffs,
    identify_equilibria
)


def test_get_coeff_mtx():
    r = [
        ({'A': 1},       {'B': 1}),
        ({'A': 1, 'B': 1}, {'C': 2})
    ]
    A = get_coeff_mtx('ABC', r)
    Aref = np.array([
        [-1, -1],
        [1, -1],
        [0,  2]
    ])
    assert np.allclose(A, Aref)


def test_decompose_yield_into_rate_coeffs():
    gamma_yields = OrderedDict([
        ('OH-', 0.5),
        ('H2O2', 0.7),
        ('OH', 2.7),
        ('H2', 0.45),
        ('H', 0.66),
        ('H+', 3.1),
        ('HO2', 0.02),
        ('e-aq', 2.6),
        # H2O: -4.64
    ])

    stoichs = [
        ({'H2O': 1}, {'H+': 1, 'OH-': 1}),
        ({'H2O': 1}, {'H+': 1, 'e-aq': 1, 'OH': 1}),
        ({'H2O': 2}, {'H': 2, 'H2O2': 1}),
        ({'H2O': 2}, {'H2': 1, 'H2O2': 1}),
        ({'H2O': 2}, {'H2': 1, 'OH': 2}),
        ({'H2O': 4}, {'H2': 3, 'HO2': 2}),
    ]

    k = decompose_yield_into_rate_coeffs(gamma_yields, stoichs)
    k_ref = [0.5, 2.6, 0.33, 0.37, 0.05, 0.01]

    assert np.allclose(k, k_ref)

    G_H2O = sum(-reac['H2O']*k[i] for
                i, (reac, prod) in enumerate(stoichs))

    assert abs(G_H2O+4.64) < 1e-3


def test_decompose_yield_into_rate_coeffs_2():
    yields = OrderedDict([
        ('B', 3.0),
        ('C', 24.0)
    ])
    stoichs = [
        ({'A': 2}, {'B': 1, 'C': 1}),
        ({'A': 1}, {'C': 3})
    ]
    k = decompose_yield_into_rate_coeffs(yields, stoichs)
    k_ref = [3, 7]

    rtol = 1e-12
    for a, b in zip(k, k_ref):
        assert abs(a-b) < abs(a*rtol)


def test_identify_equilibria():
    assert identify_equilibria([[0, 0], [1]], [[1], [0, 0]]) == set([(0, 1)])
    assert identify_equilibria([
        [0], [0, 1], [2, 3]
    ], [
        [1], [2, 3], [0, 1]
    ]) == set([(1, 2)])
