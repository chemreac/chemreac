# -*- coding: utf-8 -*-

from chemreac.util.stoich import (
    identify_equilibria
)


def test_identify_equilibria():
    assert identify_equilibria([[0, 0], [1]], [[1], [0, 0]]) == set([(0, 1)])
    assert identify_equilibria([
        [0], [0, 1], [2, 3]
    ], [
        [1], [2, 3], [0, 1]
    ]) == set([(1, 2)])
    assert identify_equilibria([
        [0, 1], [2], [3, 4]
    ], [
        [2], [3, 4], [2]
    ]) == set([(1, 2)])
