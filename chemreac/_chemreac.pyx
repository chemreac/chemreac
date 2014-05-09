# -*- coding: utf-8 -*-
# distutils: language = c++

import numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from *:
    ctypedef unsigned int uint

cdef extern from "chemreac.h" namespace "chemreac":
    cdef cppclass ReactionDiffusion:
        # (Private)
        double * D_weight

        const uint n, N, nr, nstencil
        const bool logy, logt, lrefl, rrefl
        const vector[vector[uint]] stoich_reac
        vector[vector[uint]] stoich_actv
        const vector[vector[uint]] stoich_prod
        vector[double] k
        vector[double] D
        const vector[double] x
        vector[vector[double]] bin_k_factor
        vector[uint] bin_k_factor_span
        double * xc

        ReactionDiffusion(uint,
                          const vector[vector[uint]],
                          const vector[vector[uint]],
                          vector[double],
                          uint,
                          vector[double],
                          const vector[double],
                          vector[vector[uint]],
                          vector[vector[double]],
                          vector[uint],
                          int,
                          bool,
                          bool,
                          uint,
                          bool,
                          bool) except +
        void f(double, const double * const, double * const)
        void dense_jac_rmaj(double, const double * const, double * const, int)
        void dense_jac_cmaj(double, const double * const, double * const, int)
        void banded_padded_jac_cmaj(double, const double * const, double * const, int)
        void banded_packed_jac_cmaj(double, const double * const, double * const, int)

        void per_rxn_contrib_to_fi(double, const double * const, uint, double * const)
        int get_geom_as_int()

        uint _stencil_bi_lbound(uint)
        uint _xc_bi_map(uint)


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


cdef class PyReactionDiffusion:
    cdef ReactionDiffusion *thisptr
    cdef public vector[double] k_err, D_err
    cdef public list names, tex_names

    def __cinit__(self,
                  uint n,
                  vector[vector[uint]] stoich_reac,
                  vector[vector[uint]] stoich_prod,
                  vector[double] k,
                  uint N,
                  vector[double] D,
                  vector[double] x,
                  vector[vector[uint]] stoich_actv,
                  vector[vector[double]] bin_k_factor,
                  vector[uint] bin_k_factor_span,
                  int geom,
                  bint logy,
                  bint logt,
                  uint nstencil=3,
                  bint lrefl=True,
                  bint rrefl=True,
              ):
        self.thisptr = new ReactionDiffusion(
            n, stoich_reac, stoich_prod, k, N,
            D, x, stoich_actv, bin_k_factor,
            bin_k_factor_span, geom, logy, logt, nstencil,
            lrefl, rrefl)

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
            return np.asarray(self.thisptr.D)
        def __set__(self, vector[double] D):
            assert len(D) == self.n
            self.thisptr.D = D

    property x:
        def __get__(self):
            return np.asarray(self.thisptr.x)

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

    property lrefl:
        def __get__(self):
            return self.thisptr.lrefl

    property rrefl:
        def __get__(self):
            return self.thisptr.rrefl

    # Extra convenience
    property ny:
        def __get__(self):
            return self.N*self.n

    def per_rxn_contrib_to_fi(self, double t, double[::1] y, int si, double[::1] out):
        self.thisptr.per_rxn_contrib_to_fi(t, &y[0], si, &out[0])

    property _xc:
        def __get__(self):
            return fromaddress(<long>self.thisptr.xc, (self.N+self.thisptr.nstencil-1,))

    property xcenters:
        def __get__(self):
            return fromaddress(<long>(&self.thisptr.xc[(self.thisptr.nstencil-1)//2]), (self.N,))

    # (Private)
    property D_weight:
        def __get__(self):
            return fromaddress(<long>self.thisptr.D_weight,
                                              (self.thisptr.N*self.thisptr.nstencil,))

    def _stencil_bi_lbound(self, uint bi):
        return self.thisptr._stencil_bi_lbound(bi)

    def _xc_bi_map(self, uint xci):
        return self.thisptr._xc_bi_map(xci)
