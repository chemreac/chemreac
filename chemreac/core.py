# -*- coding: utf-8 -*-

"""
core
====
Provides core functinoality: ReactionDiffusion
"""
import os
import numpy as np


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


class ReactionDiffusion(CppReactionDiffusion, ReactionDiffusionBase):
    """
    Object representing the numerical model, with callbacks for evaluating
    derivatives and jacobian.

    The instance provides methods:

    f(t, y, fout)
    dense_jac_rmaj(t, y, Jout)
    dense_jac_cmaj(t, y, Jout)
    banded_jac_cmaj(t, y, Jout)
    banded_packed_jac_cmaj(t, y, Jout)

    some of which are used by chemreac.integrate.integrate_scipy

    In addition error estimates (if provided by user) are stored as:
    - k_err
    - D_err

    Additional convenience attributes (not used by underlying C++ class):
    - substance_names
    - substance_tex_names

    Parameters
    ----------
    n: integer
        number of species
    stoich_reac: list of lists of integer indices
        reactant index lists per reaction.
    stoich_prod: list of lists of integer indices
        product index lists per reaction.
    k: 1-dimensional array
        reaction rate coefficients (if 2-tuples,
        assumed (val, stddev) pairs)
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
    bin_k_factor: sequence of sequences of floats
        per compartment modulation of rate coefficients
    bin_k_factor_span: sequence of integers
        spans over reaction indices affected by bin_k_factor
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
    eps: float
        relative permitivity of medium (dielectric constant)
    xscale: float
        use internal scaling of length (default: 1.0)
        (finite difference scheme works best for step-size ~1)
    """
    # not used by C++ class
    extra_attrs = ['k_err', 'D_err', 'substance_names', 'substance_tex_names']

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
                mobility=None, x=None, stoich_actv=None, bin_k_factor=None,
                bin_k_factor_span=None, geom=FLAT, logy=False, logt=False,
                logx=False, nstencil=None, lrefl=True, rrefl=True,
                auto_efield=False, surf_chg=(0.0, 0.0), eps=1.0,
                xscale=1.0, **kwargs):
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
            mobility = list([0]*n)
        if N > 1:
            assert n == len(D)
            assert n == len(z_chg)
            assert n == len(mobility)
        else:
            D = D or list([0]*n)

        k_val = []
        k_err = []
        D_val = []
        D_err = []
        for inp, val_lst, err_lst in [(k, k_val, k_err), (D, D_val, D_err)]:
            for entry in inp:
                try:
                    val, err = entry
                except TypeError:
                    assert isinstance(entry, float) or isinstance(entry, int)
                    val, err = entry, 0
                val_lst.append(val)
                err_lst.append(err)

        if x is None:
            x = 1.0

        if isinstance(x, float) or isinstance(x, int):
            _x = np.linspace(0, x, N+1)
        else:
            if len(x) == N+1:
                # monotonic:
                assert all([x[i+1] > x[i] for i in range(len(x)-1)])
                _x = x
            elif len(x) == 2:
                _x = np.linspace(x[0], x[1], N+1)
            else:
                raise ValueError("Don't know what to do with len(x) == %d" %
                                 len(x))

        if stoich_actv is None:
            _stoich_actv = list([[]]*len(stoich_reac))
        else:
            _stoich_actv = stoich_actv
        assert len(_stoich_actv) == len(stoich_reac)

        assert len(stoich_reac) == len(stoich_prod) == len(k)
        assert geom in (FLAT, CYLINDRICAL, SPHERICAL)

        # Handle bin_k_factor
        if bin_k_factor is None:
            if bin_k_factor_span is None:
                bin_k_factor_span = []
            bin_k_factor = []
        else:
            assert bin_k_factor_span is not None
            assert len(bin_k_factor) == N
            assert all([len(x) == len(bin_k_factor_span) for
                        x in bin_k_factor])
            assert all([x >= 0 for x in bin_k_factor_span])

        rd = super(ReactionDiffusion, cls).__new__(
            cls, n, stoich_reac, stoich_prod, k_val, N,
            D_val, z_chg, mobility, _x, _stoich_actv, bin_k_factor,
            bin_k_factor_span, geom, logy, logt, logx,
            nstencil, lrefl, rrefl, auto_efield, surf_chg, eps, xscale
        )
        rd.k_err = k_err
        rd.D_err = D_err

        for attr in cls.kwarg_attrs:
            if attr in kwargs:
                setattr(rd, '_' + attr, kwargs.pop(attr))
        if kwargs:
            raise KeyError("Unkown kwargs: ", kwargs.keys())
        return rd
