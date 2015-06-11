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
from functools import total_ordering
from itertools import chain
from operator import itemgetter
import weakref

import numpy as np
import quantities as pq

from chemreac import ReactionDiffusion
from chemreac.util.physchem import electrical_mobility_from_D

# TODO: somewhere: add check for massconservation and charge conservation

dm = pq.UnitQuantity('decimeter',  pq.m / 10.0,  symbol='dm')
molar = pq.UnitQuantity('molar',  pq.mole / dm ** 3,  symbol='M')


@total_ordering
class Substance(object):
    """
    Substance class to represent a chemical speices.

    Parameters
    ==========
    name: string
        unique string representation e.g. "H2O", "CNO-", "OCN-"
    charge: integer
        charge of substance
    mass: float
        molar mass (default None)
    formula: e.g. periodictable.formulas.Formula instance
        optional, if formula instance provides `mass` attribute it is
        used as mass in the case mass=None
    tex_name: string
        optional, TeX formated string, e.g. '$\mathrm{OH^{-}}$'
    multiplicity: integer
        optional, 1 for singlet, 2 for doublet...
    D: float (optional)
        diffusion coefficent, for now: isothermal, isotropic and only
        for one medium.  default: 0.0
    **kwargs:
        additional freely chosen attributes

    Attributes
    ==========
    all_substances
        dictionary (name, insatnce) of all Substance instances.

    Examples
    ========
    >>> Substance(name='H2O', charge=0, tex_name=r'$\mathrm{H_{2}O}$', pKa=14)
    <Substance 'H2O'>
    >>> Substance.all_substances['H2O']
    <Substance 'H2O'>
    >>> 'H2O' in Substance.all_substances
    True

    """
    # weakref => If instance is deleted GC can kill it.
    # We manually keep track och instances in order to ensure unique names
    all_substances = weakref.WeakValueDictionary()

    def __init__(self, name, charge=None, mass=None, formula=None,
                 tex_name=None, multiplicity=None, D=0.0, **kwargs):
        self.name = name
        self.charge = charge
        if mass is None and hasattr(formula, 'mass'):
            mass = formula.mass
        self.mass = mass
        self.formula = formula
        self.tex_name = tex_name
        self.multiplicity = multiplicity
        self.D = D
        self.__dict__.update(kwargs)

        if name in self.__class__.all_substances:
            colliding_occurance = self.__class__.all_substances[name]
            if not self == colliding_occurance:
                raise KeyError(
                    'Substance name already exists: ' + name + ' id=' +
                    str(id(self.__class__.all_substances[name])))
        else:
            self.__class__.all_substances[name] = self

    def __repr__(self, ):
        return "<" + self.__class__.__name__ + " '" + self.name + "'>"

    # Thanks to total_ordering it is sufficient to specify eq and lt
    def __eq__(self, other):
        """
        Equality of Substance instances is solely determined from .name
        """
        # all_substances dictionary lookup in __init__ should
        # prevent collisions
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name

    def get_mobility(self, Temp, kB=None):
        """ See ``chemreac.util.physchem.electrical_mobility_from_D`` """
        return electrical_mobility_from_D(self.D, self.charge, Temp, kB)


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
    >>> mk_sn_dict_from_names(
    ...     'ABCD', D=[0.1, 0.2, 0.3, 0.4]) # doctest: +NORMALIZE_WHITESPACE
    OrderedDict([('A', <Substance 'A'>), ('B', <Substance 'B'>),
    ('C', <Substance 'C'>), ('D', <Substance 'D'>)])

    """
    kwargs_list = []
    for i in range(len(names)):
        d = {}
        for k, v in kwargs.items():
            d[k] = v[i]
        kwargs_list.append(d)

    return OrderedDict([(s, Substance(s, **kwargs_list[i])) for i, s
                        in enumerate(names)])


class Reaction(object):
    """
    Reaction with kinetics governed by the law of mass-action.
    Example:

        A + R --> A + P; r = k*A*R

    Also supports

        5*C1 + C2 --> B; r = k*C1*C2

    by specifying active reactants C1, C2 and inactive reaktants 4*C1.

    reactants and products are dictionaries with substance names
    as keys and positive integers giving their stoichiometric coeffecients
    as values

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
    """

    # all_instances = weakref.WeakSet()

    @property
    def reactants(self):
        d = defaultdict(int)
        for k, v in chain(self.inactv_reac.items(),
                          self.active_reac.items()):
            d[k] += v
        return d

    def __init__(self, active_reac, products, inactv_reac=None,
                 k=None, T=None, Ea=None, A=None, ref=None, name=None):
        self.active_reac = defaultdict(int, active_reac)
        self.products = defaultdict(int, products)
        self.inactv_reac = defaultdict(int, inactv_reac or {})
        self.order = sum(self.active_reac.values())
        self.k = k
        self.T = T
        self.Ea = Ea
        self.A = A
        self.ref = ref
        self.name = name

    def __str__(self):
        return self.render({})

    def render(self, names, tex=False, equilibrium=False):
        if tex:
            arrow = (' $\\rightleftharpoons$ ' if equilibrium
                     else ' $\\rightarrow$ ')
        else:
            arrow = ' <-> ' if equilibrium else ' -> '
        active, inactv, prod = [[
            ((str(v)+' ') if v > 1 else '') + names.get(k, k) for
            k, v in filter(itemgetter(1), d.items())
        ] for d in (self.active_reac, self.inactv_reac, self.products)]
        fmtstr = "{}" + (" + ({})" if len(inactv) > 0 else "{}") + arrow + "{}"
        return fmtstr.format(" + ".join(active),
                             " + ".join(inactv),
                             " + ".join(prod))

    @property
    def species_names(self):
        return set(list(self.reactants.keys()) + list(self.products.keys()))

    def reactant_stoich_coeffs(self, species_names):
        return [self.reactants[n] for n in species_names]

    def product_stoich_coeffs(self, species_names):
        return [self.products[n] for n in species_names]

    # @classmethod
    # def get_reactions_with_species(cls, species_name):
    #     res = []
    #     for reaction in cls.all_instances:
    #         if (species_name in reaction.reactants.keys() or
    #            species_name in reaction.products.keys()):
    #             res.append(reaction)
    #     return res

    # def __str__(self):
    #     return ' -> '.join([
    #         ' + '.join([('' if num == 1 else str(num)) + name for
    #                     name, num in self.reactants.items() if num > 0]),
    #         ' + '.join([('' if num == 1 else str(num)) + name for
    #                     name, num in self.products.items() if num > 0])
    #     ])


class ReactionSystem(object):
    """
    Collection of reactions forming a system (model).

    Parameters
    ----------
    rxns: sequence
         sequence of :py:class:`Reaction` instances
    name: string (optional)
         Name of ReactionSystem (e.g. model name / citation key)
    substances: sequence (optional)
         Sequence of Substance instances, will be used in doing
         a sanity check and as default in method :meth:`to_ReactionDiffusion`

    Attributes
    ----------
    rxns
        sequence of :class:`Reaction` instances
    species_names
        names of occurring species
    k
        rates for rxns
    ns
        number of species
    nr
        number of reactions

    """

    def __init__(self, rxns=None, name=None, substances=None):
        self.substances = substances
        self.name = name
        self.rxns = rxns

    @property
    def rxns(self):
        return self._rxns

    @rxns.setter
    def rxns(self, val):
        self._do_sanity_check(val)
        self._rxns = val

    def _do_sanity_check(self, rxns):
        """ Check for conservation of mass and charge. """
        if self.substances is None:
            return
        for rxn in rxns:
            net_chg = 0
            net_mass = 0.0
            for reac, n in rxn.reactants.items():
                net_chg -= self.substances[reac].charge
                net_mass -= self.substances[reac].mass
            for reac, n in rxn.products.items():
                net_chg += self.substances[reac].charge
                net_mass += self.substances[reac].mass
            assert net_chg == 0
            assert abs(net_mass) < 0.01

    @classmethod
    def from_ReactionDiffusion(cls, rd):
        rxns = []
        for ri in range(rd.nr):
            rxn = rd.to_Reaction(ri)
            rxns.append(rxn)
        return cls(rxns)

    def to_ReactionDiffusion(self, substances=None, ordered_names=None,
                             **kwargs):
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
        substs = substances or self.substances

        def _kwargs_updater(key, attr):
            if attr in kwargs:
                return
            try:
                kwargs[attr] = [getattr(substs[sn], key) for sn in ord_names]
            except AttributeError:
                pass

        if substs:
            for key, attr in zip(
                    ['D', 'mobility', 'name', 'tex_name'],
                    ['D', 'mobility', 'substance_names',
                     'substance_tex_names']):
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

    def _get_repeating_indices_list(self, attr, ordered_names):
        result = []
        for rxn in self._rxns:
            l = []
            for si, sn in enumerate(ordered_names):
                l += [si]*getattr(rxn, attr, defaultdict(int))[sn]
            result.append(l)
        return result

    def stoich_active(self, ordered_names=None):
        return self._get_repeating_indices_list(
            'active_reac', ordered_names or self.ordered_names())

    def stoich_prod(self, ordered_names=None):
        return self._get_repeating_indices_list(
            'products', ordered_names or self.ordered_names())

    def stoich_inactv(self, ordered_names=None):
        return self._get_repeating_indices_list(
            'inactv_reac', ordered_names or self.ordered_names())

    @property
    def k(self):
        """ List of rate constants """
        return [rxn.k for rxn in self._rxns]

    @property
    def ns(self):
        """ Number of species """
        return len(self.species_names)

    @property
    def nr(self):
        """ Number of reactions """
        return len(self._rxns)
