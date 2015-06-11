# -*- coding: utf-8 -*-

"""
chemreac.util.stoich
--------------------

Collects stoichiometry related functions.

"""

import numpy as np


def get_coeff_mtx(substances, stoichs):
    """
    Create a net stoichiometry matrix from reactions
    described by pairs of dictionaries.

    Parameters
    ==========
    substances: sequence of keys in stoichs dict pairs
    stoichs: sequence of pairs of dicts
        pairs of reactant and product dicts mapping substance keys
        to stoichiometric coefficients (integers)

    Returns
    =======
    2 dimensional array of shape (len(substances), len(stoichs))

    """
    A = np.zeros((len(substances), len(stoichs)))
    for ri, sb in enumerate(substances):
        for ci, (reac, prod) in enumerate(stoichs):
            A[ri, ci] = prod.get(sb, 0) - reac.get(sb, 0)
    return A


def decompose_yield_into_rate_coeffs(yields, stoichs, atol=1e-10):
    """
    Decomposes (radiolytic) yields into linear combination of
    stoichiometric production reactions

    Ak = y

    A is (n_species x n_reactions) matrix, k is "rate coefficient", y is yields

    Parameters
    ==========
    yields: OrderedDict
        specie names as keys and yields as values
    stoichs: list of 2-dict tuples
         giving stoiciometry (1st is reactant, 2nd is products),
         dict keys must match those of `yields`

    Returns
    =======
    1-dimensional array of effective rate coefficients.

    """
    # Sanity check:
    for ys in yields.keys():
        present = False
        for reac, prod in stoichs:
            if ys in reac or ys in prod:
                present = True
        assert present

    sbstncs = yields.keys()
    y = np.array(list(yields.values()))
    A = get_coeff_mtx(sbstncs, stoichs)
    k, residuals, rank, s = np.linalg.lstsq(A, y)
    if len(residuals) > 0:
        if np.any(residuals > atol):
            raise ValueError("atol not satisfied")
    return k


# def get_reaction_orders(stoich_reac, stoich_actv=None):
#     """
#     Return the order of the reactions (assuming mass-action
#     behaviour).

#     Parameters
#     ----------
#     stoich_reac: list of lists of integer indices
#     stoichs: list of lists of integer indices (optional)

#     Returns
#     -------
#     iterable of integers corresponding to the total reaction orders

#     See also
#     --------
#     :class:`chemreac.core.ReactionDiffusion`

#     """
#     res = []
#     if stoich_actv is None:
#         stoich_actv = [[]]*len(stoich_reac)
#     for reac, actv in zip(stoich_reac, stoich_actv):
#         if actv == []:
#             actv = reac
#         res.append(len(actv))
#     return res


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
