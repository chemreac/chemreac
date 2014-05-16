# -*- coding: utf-8 -*-
import numpy as np

DENSE, BANDED, SPARSE = range(3)
FLAT, CYLINDRICAL, SPHERICAL = range(3)
Geom_names = {FLAT: 'Flat', CYLINDRICAL: 'Cylindrical', SPHERICAL: 'Spherical'}
GEOM_ORDER = ('Flat', 'Cylindrical', 'Spherical')

from ._chemreac import CppReactionDiffusion


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
            {self.names[i]: self.stoich_reac[ri].count(i) for\
             i in range(self.n)},
            {self.names[i]: self.stoich_prod[ri].count(i) for\
             i in range(self.n)},
            {self.names[i]: self.stoich_actv[ri].count(i) for\
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
            return np.zeros((self.n*2 + 1 + rpad, self.n*self.N + cpad), order=order)
        else:
            return np.zeros((self.n*self.N + rpad, self.n*self.N + cpad), order=order)

    @property
    def ny(self):
        return self.N*self.n

class ReactionDiffusion(CppReactionDiffusion, ReactionDiffusionBase):
    # not used by C++ class
    extra_attrs = ['k_err', 'D_err', 'names', 'tex_names']

    # subset of extra_attrs optionally passed by user
    kwarg_attrs = ['names', 'tex_names', 'xscale']

    def __new__(cls, n, stoich_reac, stoich_prod, k, N=0, D=None, x=None,
                stoich_actv=None, bin_k_factor=None, bin_k_factor_span=None,
                geom=FLAT, logy=False, logt=False, nstencil=None, lrefl=True,
                rrefl=True, **kwargs):
        """
        Arguments:
        -`n`: number of species
        -`stoich_reac`: list of reactant index lists per reaction.
        -`stoich_prod`: list of product index lists per reaction.
        -`k`: array of reaction rate coefficients (if 2-tuples, assumed (val, stddev) pairs)
        -`N`: number of compartments (default: 1 if x==None else len(x)-1)
        -`D`: diffusion coefficients (of length n)
        -`x`: compartment boundaries (of length N+1), default: linspace(1,2, N+1)
        -`stoich_actv`: list of ACTIVE reactant index lists per reaction.n
        -`bin_k_factor`: per compartment modulation of rate coefficients
        -`bin_k_factor_span`: spans over reactions affected by bin_k_factor
        -`geom`: any of (FLAT, SPHERICAL, CYLINDRICAL)
        -`logy`: f and *_jac_* routines operate on log(concentration)
        -`logt`: f and *_jac_* routines operate on log(time)
        -`nstencil`: number of points used in finite difference scheme
        -`lrefl`: reflective left boundary (default: True)
        -`rrefl`: reflective right boundary (default: True)

        Optional key-word arguments:
        -`xscale`: use internal scaling of length
                   (finite difference scheme works best for step-size ~1)

        The instance provides methods:

        f(t, y, fout)
        dense_jac_rmaj(t, y, Jout)
        dense_jac_cmaj(t, y, Jout)
        banded_jac_cmaj(t, y, Jout)
        banded_packed_jac_cmaj(t, y, Jout)

        some of which are used by chemreac.integrate.run

        In addition error estimates (if provided by user) are stored as:
        - k_err
        - D_err

        Additional convenience attributes (not used by underlying C++ class):
        - names
        - tex_names
        """

        if N == 0:
            if x == None:
                N = 1
            else:
                N = len(x)-1
        if N < nstencil:
            raise ValueError("N must be >= nstencil")

        if N > 1:
            assert n == len(D)
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

        if x == None:
            x = 1.0

        if isinstance(x, float) or isinstance(x, int):
            _x = np.linspace(1, 2, N+1)
        else:
            assert len(x) == N+1
            # monotonic:
            assert all([x[i+1]>x[i] for i in range(len(x)-1)])
            _x = x

        if stoich_actv == None:
            _stoich_actv = list([[]]*len(stoich_reac))
        else:
            _stoich_actv = stoich_actv
        assert len(_stoich_actv) == len(stoich_reac)

        assert len(stoich_reac) == len(stoich_prod) == len(k)
        assert geom in (FLAT, CYLINDRICAL, SPHERICAL)

        # Handle bin_k_factor
        if bin_k_factor == None:
            if bin_k_factor_span == None:
                bin_k_factor_span = []
            bin_k_factor = []
        else:
            assert bin_k_factor_span != None
            assert len(bin_k_factor) == N
            assert all([len(x) == len(bin_k_factor_span) for x in bin_k_factor])
            assert all([x >= 0 for x in bin_k_factor_span])

        nstencil = nstencil or (1 if N == 1 else 3)

        xscale = kwargs.pop('xscale', 1.0)

        rd = super(ReactionDiffusion, cls).__new__(
            cls, n, stoich_reac, stoich_prod, k_val, N,
            xscale**2*np.array(D_val), xscale*np.array(_x),
            _stoich_actv, bin_k_factor, bin_k_factor_span, geom, logy,
            logt, nstencil, lrefl, rrefl
        )
        rd.xscale = xscale
        rd.k_err = k_err
        rd.D_err = D_err

        for attr in cls.kwarg_attrs:
            if attr in kwargs:
                setattr(rd, attr, kwargs.pop(attr))
        if kwargs:
            raise KeyError("Unkown kwargs: ", kwargs.keys())
        return rd
