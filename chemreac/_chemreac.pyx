# -*- coding: utf-8 -*-
# distutils: language = c++

from libcpp.vector cimport vector

cdef extern from "chemreac.h" namespace "chemreac":
    cdef cppclass ReactionDiffusion:
        int n, N, nr, geom
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
                          int) except +
        void f(double, const double * const, double * const)
        void dense_jac_rmaj(double, const double * const, double * const, int)
        void dense_jac_cmaj(double, const double * const, double * const, int)
        void banded_padded_jac_cmaj(double, const double * const, double * const, int)
        void banded_packed_jac_cmaj(double, const double * const, double * const, int)


cdef class PyReactionDiffusion:
    cdef ReactionDiffusion *thisptr

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
              ):
        self.thisptr = new ReactionDiffusion(
            n, stoich_reac, stoich_prod, k, N,
            D, x, stoich_actv, bin_k_factor, bin_k_factor_span, geom)

    def __dealloc__(self):
        del self.thisptr

    def f(self, double t, double [::1] y, double [::1] fout):
        assert y.size == fout.size # OPTIMIZE AWAY
        self.thisptr.f(t, &y[0], &fout[0])

    def dense_jac_rmaj(self, double t, double [::1] y,
                       double [:, ::1] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*self.N
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.dense_jac_rmaj(t, &y[0], &Jout[0,0], Jout.shape[1])

    def dense_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*self.N
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.dense_jac_cmaj(t, &y[0], &Jout[0,0], Jout.shape[0])

    def banded_padded_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*4
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.banded_padded_jac_cmaj(t, &y[0], &Jout[0,0],
                                     Jout.shape[0])

    def banded_packed_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        assert y.size >= self.n*self.N
        assert Jout.shape[0] >= self.n*3
        assert Jout.shape[1] >= self.n*self.N
        self.thisptr.banded_packed_jac_cmaj(
            t, &y[0], &Jout[0,0], Jout.shape[0])

    property n:
        def __get__(self): return self.thisptr.n

    property N:
        def __get__(self): return self.thisptr.N

    property nr:
        def __get__(self): return self.thisptr.nr

    property stoich_reac:
        def __get__(self): return self.thisptr.stoich_reac
        def __set__(self, vector[vector[int]] stoich_reac):
            self.thisptr.stoich_reac = stoich_reac

    property stoich_prod:
        def __get__(self): return self.thisptr.stoich_prod
        def __set__(self, vector[vector[int]] stoich_prod):
            self.thisptr.stoich_prod = stoich_prod

    property stoich_actv:
        def __get__(self): return self.thisptr.stoich_actv
        def __set__(self, vector[vector[int]] stoich_actv):
            self.thisptr.stoich_actv = stoich_actv

    property k:
        def __get__(self): return self.thisptr.k
        def __set__(self, vector[double] k): self.thisptr.k = k

    property D:
        def __get__(self): return self.thisptr.D
        def __set__(self, vector[double] D): self.thisptr.D = D

    property x:
        def __get__(self): return self.thisptr.x
        def __set__(self, vector[double] x): self.thisptr.x = x

    # Extra convenience property
    property ny:
        def __get__(self): return self.N*self.n
