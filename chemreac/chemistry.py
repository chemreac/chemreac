# -*- coding: utf-8 -*-
"""
chemreac.chemistry
==================

This module collects classes useful for describing substances,
reactions and reaction systems. The classes have methods to help
with consistent low-level conversion to numerical parameters of
the model.

"""
from __future__ import print_function, division, absolute_import

from collections import defaultdict, OrderedDict
from itertools import chain
from operator import itemgetter
import weakref

import quantities as pq

from chemreac import ReactionDiffusion
from chemreac.util.physchem import electrical_mobility_from_D

from chempy.chemistry import Substance as _Substance
from chempy.chemistry import Reaction as _Reaction
from chempy.chemistry import ReactionSystem as _ReactionSystem


class Substance(_Substance):

    def get_mobility(self, Temp, **kwargs):
        """ See ``chemreac.util.physchem.electrical_mobility_from_D`` """
        return electrical_mobility_from_D(self.D, self.charge, Temp, **kwargs)


def mk_sn_dict_from_names(names, **kwargs):
    """
    Convenience function to generate a OrderedDict of Substance
    instances from a sequence of names and corresponding sequences
    of kwargs to Substance class.

    Parameters
    ----------
    names: sequence of strings
        names of substances
    **kwargs:
        sequences of corresponding keyword arguments

    Examples
    --------
    >>> d = mk_sn_dict_from_names(
    ...     'ABCD', D=[0.1, 0.2, 0.3, 0.4])
    >>> d  # doctest: +NORMALIZE_WHITESPACE
    OrderedDict([('A', <Substance 'A'>), ('B', <Substance 'B'>),
    ('C', <Substance 'C'>), ('D', <Substance 'D'>)])
    >>> d['A'].name
    'A'
    """
    kwargs_list = []
    for i in range(len(names)):
        d = {}
        other_properties = {}
        for k, v in kwargs.items():
            if k in Substance.attrs:
                d[k] = v[i]
            else:
                other_properties[k] = v[i]
        d['other_properties'] = other_properties
        kwargs_list.append(d)

    return OrderedDict([(s, Substance(s, **kwargs_list[i])) for i, s
                        in enumerate(names)])


class Reaction(_Reaction):
    """
    Reaction with kinetics governed by the law of mass-action.
    Example:

        A + R --> A + P; r = k*A*R

    Also supports

        5*C1 + C2 --> B; r = k*C1*C2

    by specifying active reactants C1, C2 and inactive reaktants 4*C1.

    reactants and products are dictionaries with substance names
    as keys and positive integers giving their stoichiometric coefficients
    as values.

    rate constant i either given as k (T optional as validity info)
    or as Ea and A for use in the Arrhenius equation

    and ref contains optional information on origin of data.

    along the same lines `name` is possible to use if the reaction
    is known under a certain name, e.g. "H2O2 auto-decomposition"

    Parameters
    ----------
    active_reac: dict
        dictionary mapping substance name (string) to stoichiometric
        coefficient (integer) of reactant, these affect rate expression.
    products: dict
        dictionary mapping substance name (string) to stoichiometric
        coefficient (integer)
    inactv_reac: dict (optional)
        Same as active_reac but does not affect rate expression.
    k: float
        rate coefficient
    T: float
        absolute temperature
    Ea: float
        activation energy
    A: float
        preexponential prefactor (Arrhenius type eq.)
    ref: string (optional)
        Reference key
    name: string (optional)
        Descriptive name of reaction

    Examples
    --------
    >>> d = mk_sn_dict_from_names('AB')
    >>> dimerization = Reaction({d['A']: 2}, {d['B']: 1})
    >>> str(dimerization)
    '2 A -> B'
    """

    def reactants(self):  # use net_stoich instead
        d = defaultdict(int)
        for k, v in chain(self.inactv_reac.items(),
                          self.active_reac.items()):
            d[k] += v
        return d

    # @property
    # def species_names(self):
    #     return set(list(self.reactants.keys()) + list(self.products.keys()))

    # def reactant_stoich_coeffs(self, species_names):
    #     return [self.reactants[n] for n in species_names]

    # def product_stoich_coeffs(self, species_names):
    #     return [self.products[n] for n in species_names]


class ReactionSystem(_ReactionSystem):

    @classmethod
    def from_ReactionDiffusion(cls, rd):
        rxns = []
        for ri in range(rd.nr):
            rxn = rd.to_Reaction(ri)
            rxns.append(rxn)
        return cls(rxns)

    def to_ReactionDiffusion(self, ordered_names=None, **kwargs):
        """
        Creates a :class:`ReactionDiffusion` instance from ``self``.

        Parameters
        ----------
        substances: sequence of Substance instances
            pass to override self.substances (optional)
        ordered_names: sequence of names
            pass to override self.ordered_names()
        \*\*kwargs:
            Keyword arguments passed on to :class:`ReactionDiffusion`
        """
        ord_names = ordered_names or self.ordered_names()

        def _kwargs_updater(key, attr):
            if attr in kwargs:
                return
            try:
                kwargs[attr] = [getattr(self.substances[sn], key) for sn in ord_names]
            except AttributeError:
                pass

        for key, attr in zip(
                ['D', 'mobility', 'name', 'latex_name'],
                ['D', 'mobility', 'substance_names',
                 'substance_latex_names']):
            _kwargs_updater(key, attr)

        return ReactionDiffusion(
            self.ns, self.stoich_active(ordered_names),
            self.stoich_prod(ordered_names), self.k,
            stoich_inactv=self.stoich_inactv(ordered_names), **kwargs)

    @property
    def species_names(self):
        return set.union(*tuple(rxn.species_names for rxn in self._rxns))

    def ordered_names(self):
        return sorted(self.species_names)

    @property
    def params(self):
        """ List of rate constants """
        return [rxn.param for rxn in self._rxns]
