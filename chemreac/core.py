# -*- coding: utf-8 -*-
"""
chemreac.core
=============
In ``chemreac.core`` you will find :py:class:`ReactionDiffusion` which
is the class describing the system of ODEs.

"""
from __future__ import (absolute_import, division, print_function)

from functools import reduce
from itertools import chain
from operator import add
import os

import numpy as np

from .chemistry import mk_sn_dict_from_names, ReactionSystem
from .util.pyutil import monotonic
from .units import unit_of, get_derived_unit, to_unitless, linspace
from .constants import get_unitless_constant

from ._chemreac import CppReactionDiffusion

DENSE, BANDED, SPARSE = range(3)
FLAT, CYLINDRICAL, SPHERICAL = range(3)
Geom_names = {FLAT: 'Flat', CYLINDRICAL: 'Cylindrical', SPHERICAL: 'Spherical'}
GEOM_ORDER = ('Flat', 'Cylindrical', 'Spherical')


class ReactionDiffusionBase(object):

    def to_ReactionSystem(self, substance_names):
        rxns = []
        for ri in range(self.nr):
            rxn = self.to_Reaction(ri, substance_names)
            rxns.append(rxn)
        return ReactionSystem(rxns, mk_sn_dict_from_names(substance_names))

    @classmethod
    def from_ReactionSystem(cls, rsys, ordered_names=None, state=None,
                            **kwargs):
        """
        Creates a :class:`ReactionDiffusion` instance from ``rsys``.

        Parameters
        ----------
        substances: sequence of Substance instances
            pass to override rsys.substances (optional)
        ordered_names: sequence of names
            pass to override rsys.ordered_names()
        state: object
            used to evaluate callable ``Reaction.params`` in ``rsys.rxns``
        \*\*kwargs:
            Keyword arguments passed on to :class:`ReactionDiffusion`
        """
        ord_names = ordered_names or rsys.substance_names()
        for rxn in rsys.rxns:
            for key in chain(rxn.reac, rxn.prod, rxn.inact_reac):
                if key not in ord_names:
                    raise ValueError("Unkown substance name: %s" % key)

        def _kwargs_updater(key, attr):
            if attr in kwargs:
                return
            try:
                kwargs[attr] = [getattr(rsys.substances[sn], key) for sn in
                                ord_names]
            except AttributeError:
                try:
                    kwargs[attr] = [rsys.substances[sn].other_properties[key]
                                    for sn in ord_names]
                except KeyError:
                    pass

        for key, attr in zip(
                ['D', 'mobility', 'name', 'latex_name'],
                ['D', 'mobility', 'substance_names',
                 'substance_latex_names']):
            _kwargs_updater(key, attr)

        return ReactionDiffusion(
            rsys.ns,
            [reduce(add, [[i]*rxn.reac.get(k, 0) for i, k
                          in enumerate(ord_names)]) for rxn in rsys.rxns],
            [reduce(add, [[i]*rxn.prod.get(k, 0) for i, k
                          in enumerate(ord_names)]) for rxn in rsys.rxns],
            [rxn.param(state) if callable(rxn.param) else rxn.param for
             rxn in rsys.rxns],
            stoich_inact=[reduce(add, [
                [i]*(0 if rxn.inact_reac is None else
                     rxn.inact_reac.get(k, 0)) for i, k in enumerate(ord_names)
            ]) for rxn in rsys.rxns],
            **kwargs)

    def to_Reaction(self, ri, substance_names=None):
        """
        Convenience method for making a Reaction instance
        for reaction index ri
        """
        from .chemistry import Reaction
        if substance_names is None:
            substance_names = self.substance_names
        return Reaction(
            {substance_names[i]: self.stoich_active[ri].count(i) for
             i in range(self.n)},
            {substance_names[i]: self.stoich_prod[ri].count(i) for
             i in range(self.n)},
            param=self.k[ri],
            inact_reac={
                substance_names[i]: self.stoich_inact[ri].count(i) for
                i in range(self.n)})

    def alloc_fout(self):
        return np.zeros(self.n*self.N)

    def alloc_jout(self, banded=None, order='C', pad=0):
        if pad is True:
            pad = self.n*self.n_jac_diags
        if pad is False:
            pad = 0
        if banded is None:
            banded = self.N > 1
        if order == 'C':
            rpad, cpad = 0, pad
        elif order == 'F':
            rpad, cpad = pad, 0
        else:
            raise ValueError("Order must be 'C' or 'F'")

        if banded:
            nr = 2*(self.n*self.n_jac_diags) + 1 + rpad
            nc = self.n*self.N + cpad
            return np.zeros((nr, nc), order=order)
        else:
            return np.zeros((self.n*self.N + rpad, self.n*self.N + cpad),
                            order=order)

    def alloc_jout_compressed(self):
        from block_diag_ilu import alloc_compressed
        return alloc_compressed(self.N, self.n, self.n_jac_diags)

    @property
    def ny(self):
        return self.N*self.n


def get_unit(unit_registry, key):
    try:
        return get_derived_unit(unit_registry, key)
    except KeyError:
        return dict(
            radyield=unit_registry['amount']/get_unit(unit_registry, 'energy'),
            field=get_unit(unit_registry, 'energy')/(get_unit(
                unit_registry, 'time') * get_unit(unit_registry, 'length')**3)
        )[key]


class ReactionDiffusion(CppReactionDiffusion, ReactionDiffusionBase):
    """
    Object representing the numerical model, with callbacks for evaluating
    derivatives and jacobian.

    The instance provides methods:

    - :meth:`f`
    - :meth:`dense_jac_rmaj`


    - ``f(t, y, fout)``
    - ``dense_jac_rmaj(t, y, jout)``
    - ``dense_jac_cmaj(t, y, jout)``
    - ``banded_jac_cmaj(t, y, jout)``
    - ``banded_packed_jac_cmaj(t, y, jout)``

    some of which are used by chemreac.integrate.integrate_scipy

    Additional convenience attributes (not used by underlying C++ class):

    - :py:attr:`substance_names`
    - :py:attr:`substance_latex_names`

    Parameters
    ----------
    n: integer
        number of species
    stoich_active: list of lists of integer indices
        reactant index lists per reaction.
    stoich_prod: list of lists of integer indices
        product index lists per reaction.
    k: 1-dimensional array
        reaction rate coefficients
    N: integer
        number of compartments (default: 1 if x==None else len(x)-1)
    D: sequence of floats
        diffusion coefficients (of length n)
    z_chg: sequence of integers
        1-dimensional array ion charges
    mobility: sequence of floats
        mobility of ions
    x: sequence of floats or pair of flats or float, optional
        compartment boundaries (of length N+1), default: linspace(0, 1, N+1)
        if x is a pair of floats it is expanded into linspace(x[0], x[1], N+1).
        if x is a float it is expanded into linspace(0, x, N+1)
    stoich_inact: list of lists of integer indices
        list of inactive reactant index lists per reaction.n, default: []
    geom: integer
        any of (FLAT, SPHERICAL, CYLINDRICAL)
    logy: bool
        f and \*_jac_\* routines operate on log(concentration)
    logt: bool
        f and \*_jac_\* routines operate on log(time)
    logx: bool
        f and \*_jac_\* routines operate on log(space)
    nstencil: integer
        number of points used in finite difference scheme
    lrefl: bool
        reflective left boundary (default: True)
    rrefl: bool
        reflective right boundary (default: True)
    auto_efield: bool
        calculate electric field from concentrations (default: False)
    surf_chg: pair of floats
        total charge of surface (defaut: (0.0, 0.0))
    eps_rel: float
        relative permitivity of medium (dielectric constant)
    g_values: sequence of sequences of floats
        per specie yields (amount / energy) per field type
    g_value_parents: sequence of integers (optional)
        indices of parents for each g_value. -1 denotes no concentration
        dependence of g_values. default: [-1]*len(g_values)
    fields: sequence of sequence of floats
        per bin field strength (energy / (volume time)) per field type.
        May be calculated as product between density and doserate
    modulated_rxns: sequence of integers
        Indicies of reactions subject to per bin modulation
    modulation: sequence of sequences of floats
        Per bin modulation vectors for each index in modulated_rxns
    ilu_limit: float
        Requirement on (average) diagonal dominance for performing ILU
        factorization.
    n_jac_diags: int
        Number of diagonals to include in Jacobian (default: 1)
    unit_registry: dict (optional)
        default: None, see ``chemreac.units.SI_base_registry`` for an
        example.

    Attributes
    ----------
    unit_registry: dict

    """
    # not used by C++ class
    extra_attrs = ['substance_names', 'substance_latex_names']

    # subset of extra_attrs optionally passed by user
    kwarg_attrs = ['substance_names', 'substance_latex_names']

    _substance_names = None
    _substance_latex_names = None

    @property
    def substance_names(self):
        return self._substance_names or list(map(str, range(self.n)))

    @substance_names.setter
    def substance_names(self, names):
        self._substance_names = names

    @property
    def substance_latex_names(self):
        return self._substance_latex_names or list(map(str, range(self.n)))

    @substance_latex_names.setter
    def substance_latex_names(self, latex_names):
        self._substance_latex_names = latex_names

    def __new__(cls, n, stoich_active, stoich_prod, k, N=0, D=None, z_chg=None,
                mobility=None, x=None, stoich_inact=None, geom=FLAT,
                logy=False, logt=False, logx=False, nstencil=None,
                lrefl=True, rrefl=True, auto_efield=False, surf_chg=None,
                eps_rel=1.0, g_values=None, g_value_parents=None,
                fields=None,
                modulated_rxns=None,
                modulation=None,
                unit_registry=None,
                faraday=None,  # deprecated
                vacuum_permittivity=None,  # deprecated
                k_unitless=None,
                ilu_limit=None,
                n_jac_diags=-1,
                **kwargs):
        if N == 0:
            if x is None:
                N = 1
            else:
                N = len(x)-1
        nstencil = nstencil or (1 if N == 1 else 3)
        if N < nstencil:
            raise ValueError("N must be >= nstencil")

        if z_chg is None:
            z_chg = list([0]*n)
        if mobility is None:
            mobility = np.zeros(n)*get_unit(unit_registry,
                                            'electrical_mobility')
        if N > 1:
            assert n == len(D)
            assert n == len(z_chg)
            assert n == len(mobility)
        else:
            if D is None:
                D = np.zeros(n)*get_unit(unit_registry, 'diffusion')

        if x is None:
            x = 1.0*get_unit(unit_registry, 'length')

        try:
            if len(x) == N+1:
                if not monotonic(x, 1, True):
                    raise ValueError("x must be strictly positive monotonic")
                _x = x
            elif len(x) == 2:
                _x = linspace(x[0], x[1], N+1)
            elif len(x) == 1:
                raise TypeError
            else:
                raise ValueError("Don't know what to do with len(x) == %d" %
                                 len(x))
        except TypeError:
            _x = linspace(0*unit_of(x), x, N+1)

        if stoich_inact is None:
            stoich_inact = list([[]]*len(stoich_active))
        assert len(stoich_inact) == len(stoich_active)
        assert len(stoich_active) == len(stoich_prod)

        if k is not None:
            assert len(stoich_active) == len(k)
        else:
            assert len(stoich_active) == len(k_unitless)
        assert geom in (FLAT, CYLINDRICAL, SPHERICAL)

        if surf_chg is None:
            surf_chg = (0.0*get_unit(unit_registry, 'charge'),
                        0.0*get_unit(unit_registry, 'charge'))

        # Handle g_values
        if g_values is None:
            g_values = []
        else:
            if fields is not None:
                assert len(g_values) == len(fields)
            for gv in g_values:
                assert len(gv) == n
        if g_value_parents is None:
            g_value_parents = [-1]*len(g_values)

        if fields is None:
            fields = [[0.0*get_unit(unit_registry, 'field')]*N]*len(g_values)
        else:
            assert len(fields) == len(g_values)
            for fld in fields:
                assert len(fld) == N

        reac_orders = map(len, stoich_active)
        if k_unitless is None:
            k_unitless = [to_unitless(kval, kunit) for kval, kunit in
                          zip(k, cls.k_units(unit_registry, reac_orders))]
        else:
            if k is not None:
                raise ValueError(
                    "When passing k_unitless you must set k to None")

        if modulated_rxns is not None:
            if len(modulated_rxns) != len(modulation):
                raise ValueError("len(modulated_rxns) != len(modulation)")
        if modulation is not None:
            if modulated_rxns is None:
                raise ValueError("modulated_rxns missing.")
            if any(len(arr) != N for arr in modulation):
                raise ValueError("An array in modulation of size != N")

        _k = np.asarray(k_unitless)
        if _k.ndim != 1:
            raise ValueError("Rates vector has inproper dimension")

        rd = super(ReactionDiffusion, cls).__new__(
            cls, n, stoich_active, stoich_prod,
            np.asarray(k_unitless),
            N,
            to_unitless(D, get_unit(unit_registry, 'diffusion')),
            z_chg,
            to_unitless(mobility, get_unit(unit_registry,
                                           'electrical_mobility')),
            to_unitless(_x, get_unit(unit_registry, 'length')),
            stoich_inact, geom, logy, logt, logx,
            [np.asarray([to_unitless(yld, yld_unit) for yld in gv]) for gv,
             yld_unit in zip(g_values, cls.g_units(unit_registry,
                                                   g_value_parents))],
            g_value_parents,
            [to_unitless(fld, get_unit(unit_registry, 'field'))
             for fld in fields],
            modulated_rxns or [],
            modulation or [],
            nstencil, lrefl, rrefl, auto_efield,
            (to_unitless(surf_chg[0], get_unit(unit_registry, 'charge')),
             to_unitless(surf_chg[1], get_unit(unit_registry, 'charge'))),
            eps_rel,
            faraday or get_unitless_constant(unit_registry,
                                             'Faraday_constant'),
            vacuum_permittivity or get_unitless_constant(
                unit_registry, 'vacuum_permittivity'),
            ilu_limit=(float(os.environ.get('CHEMREAC_ILU_LIMIT', 1000)) if
                       ilu_limit is None else ilu_limit),
            n_jac_diags=(int(os.environ.get('CHEMREAC_N_JAC_DIAGS', 1)) if
                         n_jac_diags is -1 else n_jac_diags)
        )

        rd.unit_registry = unit_registry

        for attr in cls.kwarg_attrs:
            if attr in kwargs:
                setattr(rd, '_' + attr, kwargs.pop(attr))
        if kwargs:
            raise KeyError("Unkown kwargs: ", kwargs.keys())
        return rd

    @property
    def fields(self):
        return np.asarray(self._fields)*get_unit(self.unit_registry, 'field')

    @fields.setter
    def fields(self, value):
        self._fields = to_unitless(value, get_unit(self.unit_registry,
                                                   'field'))

    @staticmethod
    def g_units(unit_registry, g_value_parents):
        """ Forms the unit of radiolytic yield

        Parameters
        ----------
        g_value_parents: iterable of integers
            parent substance indices (-1 indicates no parent)
        """
        g_units = []
        for parent in g_value_parents:
            if parent == -1:
                g_units.append(get_unit(unit_registry, 'radyield'))
            else:
                g_units.append(get_unit(unit_registry, 'radyield') /
                               get_unit(unit_registry, 'concentration'))
        return g_units

    @property
    def g_values(self):
        return [np.asarray(gv)*gu for gv, gu in zip(
            self._g_values, self.g_units(self.unit_registry,
                                         self.g_value_parents))]

    @staticmethod
    def k_units(unit_registry, reaction_orders):
        return [get_unit(unit_registry, 'concentration')**(
            1-order)/get_unit(unit_registry, 'time')
                for order in reaction_orders]

    @property
    def k(self):
        reac_orders = map(len, self.stoich_active)
        return [kv*ku for kv, ku in zip(self._k, self.k_units(
            self.unit_registry, reac_orders))]

    @k.setter
    def k(self, value):
        reac_orders = map(len, self.stoich_active)
        self._k = [to_unitless(kv, ku) for kv, ku in zip(value, self.k_units(
            self.unit_registry, reac_orders))]

    @property
    def D(self):
        return np.asarray(self._D)*get_unit(self.unit_registry, 'diffusion')

    @D.setter
    def D(self, value):
        self._D = to_unitless(value, get_unit(self.unit_registry, 'diffusion'))

    @property
    def mobility(self):
        return np.asarray(self._mobility)*get_unit(self.unit_registry,
                                                   'mobility')

    @mobility.setter
    def mobility(self, value):
        self._mobility = to_unitless(value, get_unit(self.unit_registry,
                                                     'mobility'))
