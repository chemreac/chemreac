# -*- coding: utf-8 -*-
# distutils: language = c++

import numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "chemreac.h" namespace "chemreac":
    cdef cppclass ReactionDiffusion:
        int n, N, nr, geom
        bool logy, logt
        vector[vector[int]] stoich_reac
        vector[vector[int]] stoich_actv
        vector[vector[int]] stoich_prod
        vector[double] k
        vector[double] D
        vector[double] x
        vector[vector[double]] bin_k_factor
        vector[int] bin_k_factor_span

        ReactionDiffusion(int,
                          vector[vector[int]],
                          vector[vector[int]],
                          vector[double],
                          int,
                          vector[double],
                          vector[double],
                          vector[vector[int]],
                          vector[vector[double]],
                          vector[int],
                          int,
                          bool,
                          bool) except +
        void f(double, const double * const, double * const)
        void dense_jac_rmaj(double, const double * const, double * const, int)
        void dense_jac_cmaj(double, const double * const, double * const, int)
        void banded_padded_jac_cmaj(double, const double * const, double * const, int)
        void banded_packed_jac_cmaj(double, const double * const, double * const, int)

        void per_rxn_contrib_to_fi(double, const double * const, int, double * const)
        int get_geom_as_int()


cdef class PyReactionDiffusion:
    cdef ReactionDiffusion *thisptr
    cdef public vector[double] k_err, D_err
    cdef public list names, tex_names

    def __cinit__(self,
                  int n,
                  vector[vector[int]] stoich_reac,
                  vector[vector[int]] stoich_prod,
                  vector[double] k,
                  int N,
                  vector[double] D,
                  vector[double] x,
                  vector[vector[int]] stoich_actv,
                  vector[vector[double]] bin_k_factor,
                  vector[int] bin_k_factor_span,
                  int geom,
                  bint logy,
                  bint logt,
              ):
        self.thisptr = new ReactionDiffusion(
            n, stoich_reac, stoich_prod, k, N,
            D, x, stoich_actv, bin_k_factor,
            bin_k_factor_span, geom, logy, logt)

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
        def __get__(self): return self.thisptr.n

    property N:
        def __get__(self): return self.thisptr.N

    property nr:
        def __get__(self): return self.thisptr.nr

    property geom:
        def __get__(self): return self.thisptr.get_geom_as_int()

    property stoich_reac:
        def __get__(self): return self.thisptr.stoich_reac
        # We would need to re-initialize coeff_reac, coeff_*, ...
        # def __set__(self, vector[vector[int]] stoich_reac):
        #     self.thisptr.stoich_reac = stoich_reac

    property stoich_prod:
        def __get__(self): return self.thisptr.stoich_prod
        # We would need to re-initialize coeff_reac, coeff_*, ...
        # def __set__(self, vector[vector[int]] stoich_prod):
        #     self.thisptr.stoich_prod = stoich_prod

    property stoich_actv:
        def __get__(self): return self.thisptr.stoich_actv
        # We would need to re-initialize coeff_reac, coeff_*, ...
        # def __set__(self, vector[vector[int]] stoich_actv):
        #     self.thisptr.stoich_actv = stoich_actv

    property k:
        def __get__(self): return np.asarray(self.thisptr.k)
        def __set__(self, vector[double] k):
            assert len(k) == self.nr
            self.thisptr.k = k

    property D:
        def __get__(self): return np.asarray(self.thisptr.D)
        def __set__(self, vector[double] D):
            assert len(D) == self.n
            self.thisptr.D = D

    property x:
        def __get__(self): return np.asarray(self.thisptr.x)
        def __set__(self, vector[double] x):
            assert len(x) == self.N+1
            self.thisptr.x = x

    property bin_k_factor:
        def __get__(self): return np.asarray(self.thisptr.bin_k_factor)
        def __set__(self, vector[vector[double]] bin_k_factor):
            assert len(bin_k_factor) == self.N
            self.thisptr.bin_k_factor = bin_k_factor

    property bin_k_factor_span:
        def __get__(self): return np.asarray(self.thisptr.bin_k_factor_span, dtype=np.int32)
        def __set__(self, vector[int] bin_k_factor_span):
            assert all([len(bin_k_factor_span) == len(x) for x in self.bin_k_factor])
            self.thisptr.bin_k_factor_span = bin_k_factor_span

    property logy:
        def __get__(self): return self.thisptr.logy
        def __set__(self, bool logy): self.thisptr.logy = logy

    property logt:
        def __get__(self): return self.thisptr.logt
        def __set__(self, bool logt): self.thisptr.logt = logt


    # Extra convenience
    property ny:
        def __get__(self): return self.N*self.n

    def per_rxn_contrib_to_fi(self, double t, double[::1] y, int si, double[::1] out):
        self.thisptr.per_rxn_contrib_to_fi(t, &y[0], si, &out[0])
