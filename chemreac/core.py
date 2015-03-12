# -*- coding: utf-8 -*-

"""
chemreac.core
=============
In chemreac.core you will find :py:class:`ReactionDiffusion` which
is the class describing the system of ODEs.

"""

import os
import numpy as np

from .units import unitof, get_derived_unit, to_unitless
from .util.stoich import get_reaction_orders
from .constants import get_unitless_constant

if os.environ.get('READTHEDOCS', None) == 'True':
    # On readthedocs, cannot compile extension module.
    class CppReactionDiffusion(object):
        pass  # mockup
else:
    from ._chemreac import CppReactionDiffusion, diag_data_len

DENSE, BANDED, SPARSE = range(3)
FLAT, CYLINDRICAL, SPHERICAL = range(3)
Geom_names = {FLAT: 'Flat', CYLINDRICAL: 'Cylindrical', SPHERICAL: 'Spherical'}
GEOM_ORDER = ('Flat', 'Cylindrical', 'Spherical')

# Having a Python side wrapper for our Cython Extension (CppReactionDiffusion)
# allows e.g. Jedi (Python IDE capabilities) to inspect and give help strings


class ReactionDiffusionBase(object):
    def to_Reaction(self, ri):
        """
        Convenience method for making a Reaction instance
        for reaction index ri
        """
        from .chemistry import Reaction
        return Reaction(
            {self.substance_names[i]: self.stoich_reac[ri].count(i) for
             i in range(self.n)},
            {self.substance_names[i]: self.stoich_prod[ri].count(i) for
             i in range(self.n)},
            {self.substance_names[i]: self.stoich_actv[ri].count(i) for
             i in range(self.n)},
            k=self.k[ri])

    def alloc_fout(self):
        return np.zeros(self.n*self.N)

    def alloc_jout(self, banded=True, order='C', pad=0):
        if order == 'C':
            rpad, cpad = 0, pad
        elif order == 'F':
            rpad, cpad = pad, 0
        else:
            raise ValueError("Order must be 'C' or 'F'")

        if banded:
            return np.zeros((self.n*2 + 1 + rpad, self.n*self.N + cpad),
                            order=order)
        else:
            return np.zeros((self.n*self.N + rpad, self.n*self.N + cpad),
                            order=order)

    def alloc_jout_compressed(self, ndiag):
        # TODO: ndiag from nstencil
        return np.zeros(self.n*self.n*self.N + 2*diag_data_len(
            self.N, self.n, ndiag))

    @property
    def ny(self):
        return self.N*self.n


def get_unit(units, key):
    try:
        return get_derived_unit(units, key)
    except KeyError:
        return dict(
            radyield=units['amount']/get_unit(units, 'energy'),
            field=get_unit(units, 'energy')/(get_unit(units, 'time') *
                                             get_unit(units, 'length')**3)
        )[key]


class ReactionDiffusion(CppReactionDiffusion, ReactionDiffusionBase):
    """
    Object representing the numerical model, with callbacks for evaluating
    derivatives and jacobian.

    The instance provides methods:

    - ``f(t, y, fout)``
    - ``dense_jac_rmaj(t, y, jout)``
    - ``dense_jac_cmaj(t, y, jout)``
    - ``banded_jac_cmaj(t, y, jout)``
    - ``banded_packed_jac_cmaj(t, y, jout)``

    some of which are used by chemreac.integrate.integrate_scipy

    Additional convenience attributes (not used by underlying C++ class):

    - :py:attr:`substance_names`
    - :py:attr:`substance_tex_names`

    Parameters
    ----------
    n: integer
        number of species
    stoich_reac: list of lists of integer indices
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
    stoich_actv: list of lists of integer indices
        list of ACTIVE reactant index lists per reaction.n, default: []
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
    units: dict (optional)
        default: None, see ``chemreac.units.SI_base`` for an
        example.

    Attributes
    ----------
    units: dict

    """
    # not used by C++ class
    extra_attrs = ['substance_names', 'substance_tex_names']

    # subset of extra_attrs optionally passed by user
    kwarg_attrs = ['substance_names', 'substance_tex_names']

    _substance_names = None
    _substance_tex_names = None

    @property
    def substance_names(self):
        return self._substance_names or list(map(str, range(self.n)))

    @substance_names.setter
    def substance_names(self, names):
        self._substance_names = names

    @property
    def substance_tex_names(self):
        return self._substance_tex_names or list(map(str, range(self.n)))

    @substance_tex_names.setter
    def substance_tex_names(self, tex_names):
        self._substance_tex_names = tex_names

    def __new__(cls, n, stoich_reac, stoich_prod, k, N=0, D=None, z_chg=None,
                mobility=None, x=None, stoich_actv=None, geom=FLAT,
                logy=False, logt=False, logx=False, nstencil=None,
                lrefl=True, rrefl=True, auto_efield=False, surf_chg=None,
                eps_rel=1.0, g_values=None, g_value_parents=None,
                fields=None,
                modulated_rxns=None,
                modulation=None,
                units=None,
                faraday=None,  # deprecated
                vacuum_permittivity=None,  # deprecated
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
            mobility = np.zeros(n)*get_unit(units, 'electrical_mobility')
        if N > 1:
            assert n == len(D)
            assert n == len(z_chg)
            assert n == len(mobility)
        else:
            if D is None:
                D = np.zeros(n)*get_unit(units, 'diffusion')

        if x is None:
            x = 1.0*get_unit(units, 'length')

        try:
            if len(x) == N+1:
                # monotonic:
                assert all([x[i+1] > x[i] for i in range(len(x)-1)])
                _x = x
            elif len(x) == 2:
                _x = np.linspace(x[0], x[1], N+1)
            elif len(x) == 1:
                raise TypeError
            else:
                raise ValueError("Don't know what to do with len(x) == %d" %
                                 len(x))
        except TypeError:
            _x = np.linspace(0*unitof(x), x, N+1)

        if stoich_actv is None:
            _stoich_actv = list([[]]*len(stoich_reac))
        else:
            _stoich_actv = stoich_actv
        assert len(_stoich_actv) == len(stoich_reac)

        assert len(stoich_reac) == len(stoich_prod) == len(k)
        assert geom in (FLAT, CYLINDRICAL, SPHERICAL)

        if surf_chg is None:
            surf_chg = (0.0*get_unit(units, 'charge'),
                        0.0*get_unit(units, 'charge'))

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
            fields = [[0.0*get_unit(units, 'field')]*N]*len(g_values)
        else:
            assert len(fields) == len(g_values)
            for fld in fields:
                assert len(fld) == N

        k_unitless = []
        for order, kval in zip(get_reaction_orders(stoich_reac, stoich_actv),
                               k):
            k_unitless.append(to_unitless(kval, get_unit(
                units, 'concentration')**(1-order)/get_unit(units, 'time')))
        g_units = []
        for parent in g_value_parents:
            if parent == -1:
                g_units.append(get_unit(units, 'radyield'))
            else:
                g_units.append(get_unit(units, 'radyield') /
                               get_unit(units, 'concentration'))

        rd = super(ReactionDiffusion, cls).__new__(
            cls, n, stoich_reac, stoich_prod,
            np.asarray(k_unitless),
            N,
            to_unitless(D, get_unit(units, 'diffusion')),
            z_chg,
            to_unitless(mobility, get_unit(units, 'electrical_mobility')),
            to_unitless(_x, get_unit(units, 'length')),
            _stoich_actv, geom, logy, logt, logx,
            [np.asarray([to_unitless(yld, yld_unit) for yld in gv]) for
             gv, yld_unit in zip(g_values, g_units)],
            g_value_parents,
            [to_unitless(fld, get_unit(units, 'field')) for fld in fields],
            modulated_rxns or [],
            modulation or [],
            nstencil, lrefl, rrefl, auto_efield,
            (to_unitless(surf_chg[0], get_unit(units, 'charge')),
             to_unitless(surf_chg[1], get_unit(units, 'charge'))),
            eps_rel,
            faraday or get_unitless_constant(units, 'faraday'),
            vacuum_permittivity or get_unitless_constant(
                units, 'vacuum_permittivity'),
        )

        rd.units = units

        for attr in cls.kwarg_attrs:
            if attr in kwargs:
                setattr(rd, '_' + attr, kwargs.pop(attr))
        if kwargs:
            raise KeyError("Unkown kwargs: ", kwargs.keys())
        return rd

    @property
    def fields(self):
        return np.asarray(self._fields)*get_unit(self.units, 'field')

    @fields.setter
    def fields(self, value):
        self._fields = value/get_unit(self.units, 'field')
