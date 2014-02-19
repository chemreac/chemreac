# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import

import weakref
from collections import defaultdict, OrderedDict
from functools import total_ordering

import numpy as np
import quantities as pq

from chemreac.helpers import InstanceReferenceStore
from chemreac import ReactionDiffusion

# TODO: somewhere: add check for massconservation and charge conservation

dm =  pq.UnitQuantity('decimeter',  pq.m / 10.0,  symbol = 'dm')
molar =  pq.UnitQuantity('molar',  pq.mole / dm ** 3,  symbol = 'M')


@total_ordering
class Substance(object):
    """
    sorted on name
    """
    # weakref => If instance is deleted GC can kill it.
    # Instead of subclassing InstanceReferenceStore we manually keep
    # track och instances in order to ensure unique names
    all_substances = weakref.WeakValueDictionary()

    def __init__(self,
                 name,       # Must be unique
                 charge   = None,
                 mass     = None,
                 formula  = None,
                 tex_name = None,
                 pKa      = None,
                 multiplicity = None,
                 D = None,
                 ):

        self.name = name
        self.charge = charge
        if mass == None and hasattr(formula, 'mass'):
            mass = formula.mass
        self.mass = mass
        self.formula = formula
        self.tex_name = tex_name
        self.pKa = pKa
        self.multiplicity = multiplicity
        self.D = D

        if name in self.__class__.all_substances:
            colliding_occurance = self.__class__.all_substances[name]
            if not self == colliding_occurance:
                raise KeyError('Substance name already exists: '+\
                               name+' id='+str(id(\
                                   self.__class__.all_substances[\
                                       name])))
        else:
            self.__class__.all_substances[name] = self


    def __str__(self):
        return self.name

    def __repr__(self):
        fmtstr =  "Substance('{}', {}, {}, {}, '{}', {}, {}, {})"
        return fmtstr.format(self.name, self.charge,
                             self.mass, repr(self.formula),
                             self.tex_name,
                             self.pKa, self.multiplicity, self.D)

    def __hash__(self):
        return hash(repr(self))

    # Thanks to total_ordering it is sufficient to specify eq and lt
    def __eq__(self, other):
        """
        Equality of Substance instances is solely determined from .name
        """
        # all_substances dictionary lookup in __init__ should
        # prevent collisions
        return self.name == other.name


    def __lt__(self, other):
        return self.name <  other.name




def mk_sn_dict_from_names(names, **kwargs):
    kwargs_list = [{k: v[i]} for k,v in kwargs.items() for i in range(len(names))]
    return OrderedDict([(s, Substance(s, **kwargs_list[i])) for \
                        i,s in enumerate(names)])


class Reaction(InstanceReferenceStore):
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
    or as E_a and A for use in the Arrhenius equation

    and ref contains optional information on origin of data.

    along the same lines `name` is possible to use if the reaction
    is known under a certain name, e.g. "H2O2 auto-decomposition"
    """
    def __init__(self,
                 reactants = defaultdict(int),
                 products  = defaultdict(int),
                 active_reac = None,
                 inactive_reac = None,
                 k         = None,
                 T         = None,
                 E_a       = None,
                 A         = None,
                 ref       = None,
                 name      = None,
                 ):
        # Initialize the parent:
        super(self.__class__, self).__init__()

        self._reactants = defaultdict(int)
        self._reactants.update(reactants)
        self._products  = defaultdict(int)
        self._products.update(products)
        self._active_reac = defaultdict(int)
        if active_reac:
            assert inactive_reac == None
            self._active_reac.update(active_reac)
        else:
            if inactive_reac:
                assert active_reac == None
                self._active_reac.update(reactants)
                for k, v in inactive_reac.items:
                    self._active_reac[k] -= v

        self._k     = k
        self._T     = T
        self._E_a   = E_a
        self._A     = A
        self._ref   = ref
        self._name  = name


    def __str__(self):
        return self._name if self._name else self.__repr__()


    def __repr__(self):
        attrs = ['active_reac', 'k', 'T', 'E_a', 'A', 'ref', 'name']
        pres_attrs = []
        for attr in attrs:
            if getattr(self, '_' + attr) != None:
                pres_attrs.append(attr)

        fmtcore = ', '.join(['{}'] * (2 + len(pres_attrs)))
        fmtstr =  'Reaction(' +fmtcore+ ')'
        fmtargs = [self._reactants, self._products] + \
            [attr +' = '+str(getattr(self, '_' + attr)) for attr in pres_attrs]
        return fmtstr.format(*fmtargs)


    @property
    def species_names(self):
        return set(self._reactants.keys() + self._products.keys() + \
                   (self._active_reac or {}).keys())


    def reactant_stoich_coeffs(self, species_names):
        return [self._reactants[n] for n in species_names]


    def product_stoich_coeffs(self, species_names):
        return [self._products[n] for n in species_names]


    @classmethod
    def get_reactions_with_species(cls, species_name):
        res = []
        for reaction in cls.get_all_instances():
            if species_name in reaction._reactants.keys() or \
                   species_name in reaction._products.keys():
                res.append(reaction)
        return res


class ReactionSystem(object):

    def __init__(self,
                 rxns = None,   # Reactions in system
                 T    = None, # Kelvin
                 name = None,
                 ):

        self._rxns = rxns
        self._T = T
        self._name = name

    def to_ReactionDiffusion(self, substances=None, ordered_names=None, **kwargs):
        ordered_names = ordered_names or self.ordered_names()
        if substances:
            if not 'D' in kwargs:
                kwargs['D'] = [substances[sn].D for sn in ordered_names]
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
        return self._get_repeating_indices_list('_reactants', ordered_names or self.ordered_names())


    def stoich_prod(self, ordered_names=None):
        return self._get_repeating_indices_list('_products', ordered_names or self.ordered_names())


    def stoich_actv(self, ordered_names=None):
        return self._get_repeating_indices_list('_active_reac', ordered_names or self.ordered_names())


    @property
    def k(self):
        return [rxn._k for rxn in self._rxns]


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

    def __init__(self, k_H0, derivative, ref = None):
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
        return self._k_H0 * np.exp(self._derivative \
                                   * (1 / T - 1 / self.T0))

    def get_c_at_T_and_P(self, T, P):
        return P * self.get_k_H_at_T(T)

    def get_P_at_T_and_c(self, T, c):
        return c / self.get_k_H_at_T(T)

    def plot(self, T_start, T_end):
        T = np.linspace(T_start, T_end, 100)
        k_H = self.k_H(T)
        plt.plot(T, k_H)


    def __repr__(self):
        return self.__class__.__name__ + '({}, {}, {})'.format(
            self._k_H0, self._derivative, self._ref)
