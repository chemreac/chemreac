# -*- coding: utf-8 -*-

"""
chemreac.util.stoich
--------------------

Collects stoichiometry related functions.

"""


def identify_equilibria(stoich_reac, stoich_prod):
    """
    Identify equilibria from stoichiometry

    Parameters
    ----------
    stoich_reac: iterable of iterables of integers
        per reaction iterables of specie indices for reactants
    stoich_prod: iterable of iterables of integers
        per reaction iterables of specie indices for products

    Returns
    -------
    Set of tuples of reaction indices forming equilibria

    Examples
    --------
    >>> identify_equilibria([[0,0], [1]], [[1], [0,0]]) == set([(0, 1)])
    True

    """
    equilibria = set()
    rxns = tuple(zip(stoich_reac, stoich_prod))
    for ri, (cur_reac, cur_prod) in enumerate(rxns):
        for oi, (other_reac, other_prod) in enumerate(rxns[ri+1:], start=ri+1):
            if cur_reac == other_prod and cur_prod == other_reac:
                equilibria.add((ri, oi))
    return equilibria
