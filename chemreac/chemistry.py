# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

from collections import defaultdict, OrderedDict
from functools import total_ordering
from operator import itemgetter
import weakref

import numpy as np
import quantities as pq

from chemreac import ReactionDiffusion

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
         unique string representation e.g. "CNO-", "ONC-"
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
    D: diffusion coefficent
         optional, default 0.0, for now: isothermal, isotropic and only
         for one medium.
    **kwargs:
         additional freely chosen attributes

    Examples
    ========
    >>> water = Substance(name = 'H2O',  charge = 0, formula = formula_H2O,
        tex_name = r'$\mathrm{H_{2}O}$', pKa = 14)
    <Substance 'H2O'>
    >>> Substance.all_substances['H2O']
    <Substance 'H2O'>
    >>> del water
    >>> 'H2O' in Substance.all_substances
    False


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
        return "<" + __class__.__name__ + " '" + self.name + "'>"

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


def mk_sn_dict_from_names(names, **kwargs):
    """
    Convenience function to generate a OrderedDict of Substance
    instances from a sequence of names and corresponding sequences
    of kwargs to Substance class.

    Examples
    =======
    >>> sbstncs = mk_sn_dict_from_names('ABCD', D=[0.1, 0.2, 0.3, 0.4])
    OrderedDict([('A', <Substance 'A'>), 'B', <Substance 'B'>), 'C',
    <Substance 'C'>), 'D', <Substance 'D'>)])
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

    Must honour:
      A + R --> A + P
    That is: law of massaction depend on [A]

    Also supports

     5*C1 + C2 --> B;  r=k*C1*C2

    by specifying active reactants

    reactants and products are dictionaries with substance names
    as keys and positive integers giving their stoichiometric coeffecients
    as values

    rate constant i either given as k (T optional as validity info)
    or as Ea and A for use in the Arrhenius equation

    and ref contains optional information on origin of data.

    along the same lines `name` is possible to use if the reaction
    is known under a certain name, e.g. "H2O2 auto-decomposition"
    """

    all_instances = weakref.WeakSet()

    def __init__(self, reactants, products, active_reac=None,
                 inactive_reac=None, k=None, T=None, Ea=None, A=None,
                 ref=None, name=None):
        self.all_instances.add(self)
        self.reactants = defaultdict(int)
        self.reactants.update(reactants)
        self.products = defaultdict(int)
        self.products.update(products)
        self.active_reac = defaultdict(int)
        if active_reac:
            assert inactive_reac is None
            self.active_reac.update(active_reac)
        else:
            if inactive_reac:
                assert active_reac is None
                self.active_reac.update(reactants)
                for key, val in inactive_reac.items():
                    self.active_reac[key] -= val

        self.k = k
        self.T = T
        self.Ea = Ea
        self.A = A
        self.ref = ref
        self.name = name

    def __str__(self):
        return self.render({})

    def render(self, tex_names):
        if len(tex_names) > 0:
            arrow = ' $\\rightarrow$ '
        else:
            arrow = ' -> '

        reac, prod = [[
            ((str(v)+' ') if v > 1 else '') + tex_names.get(k, k) for
            k, v in filter(itemgetter(1), d.items())
        ] for d in (self.reactants, self.products)]
        return " + ".join(reac) + arrow + " + ".join(prod)

    @property
    def species_names(self):
        return set(self.reactants.keys() + self.products.keys() +
                   (self.active_reac or {}).keys())

    def reactant_stoich_coeffs(self, species_names):
        return [self.reactants[n] for n in species_names]

    def product_stoich_coeffs(self, species_names):
        return [self.products[n] for n in species_names]

    @classmethod
    def get_reactions_with_species(cls, species_name):
        res = []
        for reaction in self.all_instances:
            if (species_name in reaction.reactants.keys() or
               species_name in reaction.products.keys()):
                res.append(reaction)
        return res


class ReactionSystem(object):
    """
    Collection of reactions forming a system (model).

    Parameters
    ==========
    rxns: sequence
         Reaction instances in ReactionSystem
    name: string
         Name of ReactionSystem (e.g. model name / citation key)
    substances: sequence
         Sequence of Substance instances, will be used in doing
         a sanity check and as default in method `to_ReactionDiffusion`
    """

    def __init__(self, rxns=None, name=None, substances=None):
        self._rxns = rxns
        self.name = name
        self.substances = substances

        self._do_sanity_check()

    @property
    def rxns(self):
        return self._rxns

    @rxns.setter
    def rxns(self, val):
        self._rxns = val
        self._do_sanity_check()

    def _do_sanity_check(self):
        if self.substances is None:
            return
        for rxn in self._rxns:
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

    def to_ReactionDiffusion(self, substances=None, ordered_names=None,
                             **kwargs):
        ord_names = ordered_names or self.ordered_names()
        substs = substances or self.substances

        def _kwargs_updater(key, attr):
            if attr in kwargs:
                return
            kwargs[attr] = [getattr(substs[sn], key) for sn in ord_names]

        if substs:
            for key, attr in zip(
                    ['D', 'name', 'tex_name'],
                    ['D', 'substance_names', 'substance_tex_names']):
                _kwargs_updater(key, attr)

        assert 'stoich_actv' not in kwargs
        return ReactionDiffusion(
            self.ns, self.stoich_reac(ordered_names),
            self.stoich_prod(ordered_names), self.k,
            stoich_actv=self.stoich_actv(ordered_names), **kwargs)

    @property
    def species_names(self):
        return set.union(*(rxn.species_names for rxn in self._rxns))

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

    def stoich_reac(self, ordered_names=None):
        return self._get_repeating_indices_list(
            'reactants', ordered_names or self.ordered_names())

    def stoich_prod(self, ordered_names=None):
        return self._get_repeating_indices_list(
            'products', ordered_names or self.ordered_names())

    def stoich_actv(self, ordered_names=None):
        return self._get_repeating_indices_list(
            'active_reac', ordered_names or self.ordered_names())

    @property
    def k(self):
        return [rxn.k for rxn in self._rxns]

    @property
    def ns(self):
        return len(self.species_names)

    @property
    def nr(self):
        return len(self._rxns)


class Henry(object):
    """
    Henry's gas constant
    """

    T0 = 298.15 * pq.kelvin

    def __init__(self, k_H0, derivative, ref=None):
        """

        Arguments:
        - `k_H0`:        Henry's constant [M/atm]
        - `derivative`: -dln(k_H)/d(1/T) [K]
        - `ref`:         Note about origin of parameters
        """
        self._k_H0 = k_H0
        self._derivative = derivative
        self._ref = ref

    def get_k_H_at_T(self, T):
        return self._k_H0 * np.exp(
            self._derivative*(1/T - 1/self.T0))

    def get_c_at_T_and_P(self, T, P):
        return P * self.get_k_H_at_T(T)

    def get_P_at_T_and_c(self, T, c):
        return c / self.get_k_H_at_T(T)
