# -*- coding: utf-8 -*-
# distutils: language = c++

import numpy as np
cimport numpy as cnp

from chemreac cimport ReactionDiffusion
from chemreac_sundials cimport direct

from libcpp.vector cimport vector

cdef extern from *:
    ctypedef unsigned int uint

DEF FLAT=0
DEF CYLINDRICAL=1
DEF SPHERICAL=2



cdef class ArrayWrapper(object):
    cdef public dict __array_interface__
    def __init__(self, **kwargs):
        self.__array_interface__ = kwargs


cdef fromaddress(address, shape, dtype=np.float64, strides=None, ro=True):
    dtype = np.dtype(dtype)
    return np.asarray(ArrayWrapper(
        data = (address, ro),
        typestr = dtype.str,
        descr = dtype.descr,
        shape = shape,
        strides = strides,
        version = 3,
    ))


cdef class CppReactionDiffusion:
    """
    Wrapper around C++ class ReactionDiffusion,
    In addition of being a thing wrapper it abstracts:
        -`xscale`: scaling of length (affects D and mobility)

    """
    cdef ReactionDiffusion *thisptr
    cdef public vector[double] k_err, D_err
    cdef public list names, tex_names
    cdef readonly double xscale

    def __cinit__(self,
                  uint n,
                  vector[vector[uint]] stoich_reac,
                  vector[vector[uint]] stoich_prod,
                  vector[double] k,
                  uint N,
                  vector[double] D,
                  vector[int] z_chg,
                  vector[double] mobility,
                  vector[double] x,
                  vector[vector[uint]] stoich_actv,
                  vector[vector[double]] bin_k_factor,
                  vector[uint] bin_k_factor_span,
                  int geom,
                  bint logy,
                  bint logt,
                  bint logx,
                  uint nstencil=3,
                  bint lrefl=True,
                  bint rrefl=True,
                  bint auto_efield=False,
                  double surf_chg=0.0,
                  double eps=1.0,
                  double xscale = 1.0,
              ):
        cdef size_t i
        for i in range(x.size()):
            x[i] *= xscale
        for i in range(D.size()):
            D[i] *= xscale**2
        self.xscale = xscale
        self.thisptr = new ReactionDiffusion(
            n, stoich_reac, stoich_prod, k, N,
            D, z_chg, mobility, x, stoich_actv, bin_k_factor,
            bin_k_factor_span, geom, logy, logt, logx, nstencil,
            lrefl, rrefl, auto_efield, surf_chg, eps)

    def __dealloc__(self):
        del self.thisptr

    def f(self, double t, double [::1] y, double [::1] fout):
        assert y.size == fout.size
        assert y.size >= self.n
        self.thisptr.f(t, &y[0], &fout[0])

    def dense_jac_rmaj(self, double t, double [::1] y,
                       double [:, ::1] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*self.N
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.dense_jac_rmaj(
            t, &y[0], &Jout[0,0], Jout.shape[1])

    def dense_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*self.N
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.dense_jac_cmaj(
            t, &y[0], &Jout[0,0], Jout.shape[0])

    def banded_padded_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*3+1
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.banded_padded_jac_cmaj(
            t, &y[0], &Jout[0,0], Jout.shape[0])

    def banded_packed_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*2+1
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.banded_packed_jac_cmaj(
            t, &y[0], &Jout[0,0], Jout.shape[0])

    def calc_efield(self, double[::1] linC):
        self.thisptr.calc_efield(&linC[0])
        return self.efield  # convenience

    def integrated_conc(self, linC):
        """
        Integrates the concentration over the volume of the system.
        Pass linear concentration "linC"
        """
        if linC.shape != (self.N,):
            raise ValueError("linC must be of length N")
        if self.geom == FLAT:
            return np.sum(np.diff(self.lin_x)*linC)
        elif self.geom == CYLINDRICAL:
            return np.sum(np.pi*np.diff(self.lin_x**2)*linC)
        elif self.geom == SPHERICAL:
            return np.sum(4*np.pi/3*np.diff(self.lin_x**3)*linC)

    property lin_x:
        def __get__(self):
            if self.logx:
                return np.exp(self.x)
            else:
                return self.x

    property n:
        def __get__(self):
            return self.thisptr.n

    property N:
        def __get__(self):
            return self.thisptr.N

    property nr:
        def __get__(self):
            return self.thisptr.nr

    property geom:
        def __get__(self):
            return self.thisptr.get_geom_as_int()

    property stoich_reac:
        def __get__(self):
            return self.thisptr.stoich_reac

    property stoich_prod:
        def __get__(self):
            return self.thisptr.stoich_prod

    property stoich_actv:
        def __get__(self):
            return self.thisptr.stoich_actv

    property k:
        def __get__(self):
            return np.asarray(self.thisptr.k)
        def __set__(self, vector[double] k):
            assert len(k) == self.nr
            self.thisptr.k = k

    property D:
        def __get__(self):
            return np.asarray(self.thisptr.D)/self.xscale**2

        def __set__(self, vector[double] D):
            cdef size_t i
            for i in range(D.size()):
                D[i] *= self.xscale**2
            assert len(D) == self.n
            self.thisptr.D = D

    property z_chg:
        def __get__(self):
            return np.asarray(self.thisptr.z_chg, dtype=np.int32)

        def __set__(self, vector[int] z_chg):
            assert len(z_chg) == self.n
            self.thisptr.z_chg = z_chg

    property mobility:
        def __get__(self):
            return np.asarray(self.thisptr.mobility)/self.xscale**2

        def __set__(self, vector[double] mobility):
            cdef size_t i
            for i in range(mobility.size()):
                mobility[i] *= self.xscale**2
            assert len(mobility) == self.n
            self.thisptr.mobility = mobility

    property x:
        def __get__(self):
            return np.asarray(self.thisptr.x)/self.xscale

    property bin_k_factor:
        def __get__(self):
            return np.asarray(self.thisptr.bin_k_factor)
        def __set__(self, vector[vector[double]] bin_k_factor):
            assert len(bin_k_factor) == self.N
            self.thisptr.bin_k_factor = bin_k_factor

    property bin_k_factor_span:
        def __get__(self):
            return np.asarray(self.thisptr.bin_k_factor_span, dtype=np.int32)
        def __set__(self, vector[uint] bin_k_factor_span):
            assert all([len(bin_k_factor_span) == len(x) for x in self.bin_k_factor])
            self.thisptr.bin_k_factor_span = bin_k_factor_span

    property logy:
        def __get__(self):
            return self.thisptr.logy

    property logt:
        def __get__(self):
            return self.thisptr.logt

    property logx:
        def __get__(self):
            return self.thisptr.logx

    property nstencil:
        def __get__(self):
            return self.thisptr.nstencil

    property lrefl:
        def __get__(self):
            return self.thisptr.lrefl

    property rrefl:
        def __get__(self):
            return self.thisptr.rrefl

    property auto_efield:
        def __get__(self):
            return self.thisptr.auto_efield

    property surf_chg:
        def __get__(self):
            return self.thisptr.surf_chg

    property eps:
        def __get__(self):
            return self.thisptr.eps

    # Extra convenience
    def per_rxn_contrib_to_fi(self, double t, double[::1] y, int si, double[::1] out):
        self.thisptr.per_rxn_contrib_to_fi(t, &y[0], si, &out[0])

    property xcenters:
        def __get__(self):
            return 1/self.xscale*fromaddress(
                <long>(&self.thisptr.xc[(self.thisptr.nstencil-1)//2]), (self.N,))

    # For debugging
    property _xc:
        def __get__(self):
            return fromaddress(<long>self.thisptr.xc, (self.N+self.thisptr.nstencil-1,))

    property D_weight:  # (Private)
        def __get__(self):
            return fromaddress(<long>self.thisptr.D_weight,
                               (self.thisptr.N*self.thisptr.nstencil,))

    def _stencil_bi_lbound(self, uint bi):
        return self.thisptr._stencil_bi_lbound(bi)

    def _xc_bi_map(self, uint xci):
        return self.thisptr._xc_bi_map(xci)

    property efield:
        def __get__(self):
            return fromaddress(<long>self.thisptr.efield, (self.thisptr.N,))
        def __set__(self, double[:] efield):
            cdef int i
            assert efield.size == self.thisptr.N
            for i in range(self.thisptr.N):
                self.thisptr.efield[i] = efield[i]


def sundials_direct(
        CppReactionDiffusion rd, vector[double] atol, double rtol,
        basestring lmm, double[::1] y0, double[::1] tout):
    cdef cnp.ndarray[cnp.float64_t, ndim=1] yout = np.empty((tout.size+1)*rd.n*rd.N)
    assert y0.size == rd.n
    print("about to call direct[]() and right now rd.auto_efield=",rd.auto_efield)
    direct[double, ReactionDiffusion](rd.thisptr, atol, rtol,
                                      {'adams': 1, 'bdf': 2}[lmm.lower()], &y0[0],
                                      tout.size, &tout[0], &yout[0])
    return yout.reshape((tout.size, rd.N, rd.n))
