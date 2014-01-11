# -*- coding: utf-8 -*-
# distutils: language = c++

from libcpp.vector cimport vector

cdef extern from "cpp_chem.hpp":
    cdef cppclass ReactionDiffusion:
        int n, N, nr, mode, geom
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
                          int) except +
        void f(double, const double * const, double * const)
        void dense_jac_rmaj(double, const double * const, double * const, int)
        void dense_jac_cmaj(double, const double * const, double * const, int)
        void banded_jac_cmaj(double, const double * const, double * const, int)
        void banded_packed_jac_cmaj(double, const double * const, double * const, int)


DEF DENSE=0
DEF BANDED=1
DEF SPARSE=2

DEF FLAT=0
DEF SPHERICAL=1
DEF CYLINDRICAL=2

cdef class PyReactionDiffusion:
    cdef ReactionDiffusion *thisptr

    def __cinit__(self,
                  int n,
                  vector[vector[int]] stoich_reac,
                  vector[vector[int]] stoich_prod,
                  vector[double] k,
                  int N = 0,
                  D = None,
                  x = None,
                  stoich_actv = None,
                  bin_k_factor = None,
                  bin_k_factor_span = None,
                  int mode=DENSE,
                  int geom=FLAT,
              ):
        cdef vector[vector[int]] _stoich_actv
        cdef vector[double] _x
        cdef vector[double] _D

        if N == 0:
            if x == None:
                N = 1
            else:
                N = len(x)-1

        if N > 1:
            assert n == len(D)
            _D = D
        else:
            _D = list([0]*n)

        x = x or 1

        if isinstance(x, float) or isinstance(x, int):
            _x = [x/float(N)*i for i in range(N+1)]
        else:
            assert len(x) == N+1
            _x = x

        if stoich_actv == None:
            _stoich_actv = list([[]]*len(stoich_reac))
        else:
            _stoich_actv = stoich_actv
        assert len(_stoich_actv) == len(stoich_reac)

        assert len(stoich_reac) == len(stoich_prod) == len(k)
        assert mode in (DENSE, BANDED, SPARSE)
        assert geom in (FLAT, SPHERICAL, CYLINDRICAL)

        # Handle bin_k_factor
        if bin_k_factor == None:
            assert bin_k_factor_span == None
            bin_k_factor = []
            bin_k_factor_span = []
        else:
            assert bin_k_factor_span != None
            assert len(bin_k_factor) == N
            assert all([len(x) == len(bin_k_factor_span) for x in bin_k_factor])
            assert all([x >= 0 for x in bin_k_factor_span])

        self.thisptr = new ReactionDiffusion(
            n, stoich_reac, stoich_prod, k, N,
            _D, _x, _stoich_actv, bin_k_factor, bin_k_factor_span, mode, geom)

    def __dealloc__(self):
        del self.thisptr

    def f(self, double t, double [::1] y, double [::1] fout):
        assert y.size == fout.size # OPTIMIZE AWAY
        self.thisptr.f(t, &y[0], &fout[0])

    def dense_jac_rmaj(self, double t, double [::1] y,
                       double [:, ::1] Jout):
        self.thisptr.dense_jac_rmaj(t, &y[0], &Jout[0,0], Jout.shape[1])

    def dense_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        self.thisptr.dense_jac_cmaj(t, &y[0], &Jout[0,0], Jout.shape[1])

    def banded_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
        self.thisptr.banded_jac_cmaj(t, &y[0], &Jout[0,0],
                                     Jout.shape[0])

    def banded_packed_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] Jout):
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
