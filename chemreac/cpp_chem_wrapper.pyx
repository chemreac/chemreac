# distutils: language = c++

from libcpp.vector cimport vector

cdef extern from "cpp_chem.hpp":
    cdef cppclass ReactionDiffusion:
        int n, N, nr
        vector[vector[int]] stoich_reac
        vector[vector[int]] stoich_prod
        vector[double] k
        vector[double] D
        int nfeval
        int njeval

        ReactionDiffusion(int, int, vector[vector[int]], vector[vector[int]],
                          vector[double], vector[double], int) except +
        void f(double, double*, double*)
        void dense_jac_rmaj(double, double*, double*, int, double, int)
        void dense_jac_cmaj(double, double*, double*, int, double, int)
        void banded_jac_cmaj(double, double*, double*, int, double, int)
        void banded_packed_jac_cmaj(double, double*, double*, int, double, int)

DEF DENSE=0
DEF BANDED=1
DEF SPARSE=2

cdef class PyReactionDiffusion:
    cdef ReactionDiffusion *thisptr

    def __cinit__(self, int n, int N,
                  vector[vector[int]] stoich_reac,
                  vector[vector[int]] stoich_prod,
                  vector[double] k,
                  vector[double] D,
                  int mode):
        assert len(stoich_reac) == len(stoich_prod) == len(k)
        assert n == len(D)
        assert mode in range(3)
        self.thisptr = new ReactionDiffusion(
            n, N, stoich_reac, stoich_prod, k, D, mode)

    def __dealloc__(self):
        del self.thisptr

    def f(self, double t, double [::1] y, double [::1] dydt):
        self.thisptr.f(t, &y[0], &dydt[0])

    def dense_jac_rmaj(self, double t, double [::1] y,
                       double [:, ::1] J, double factor, bint sub_one):
        self.thisptr.dense_jac_rmaj(
            t, &y[0], &J[0,0], J.shape[1], factor, sub_one)

    def dense_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] J, double factor, bint sub_one):
        self.thisptr.dense_jac_cmaj(
            t, &y[0], &J[0,0], J.shape[1], factor, sub_one)

    def banded_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] J, double factor, bint sub_one):
        self.thisptr.banded_jac_cmaj(
            t, &y[0], &J[0,0], J.shape[0], factor, sub_one)

    def banded_packed_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] J, double factor, bint sub_one):
        self.thisptr.banded_packed_jac_cmaj(
            t, &y[0], &J[0,0], J.shape[0], factor, sub_one)

    property n:
        def __get__(self): return self.thisptr.n
#        def __set__(self, int n): self.thisptr.n = n

    property N:
        def __get__(self): return self.thisptr.N
        def __set__(self, int N): self.thisptr.N = N

    property nr:
        def __get__(self): return self.thisptr.nr
#        def __set__(self, int nr): self.thisptr.nr = nr

    property nfeval:
        def __get__(self): return self.thisptr.nfeval
        def __set__(self, int nfeval): self.thisptr.nfeval = nfeval

    property njeval:
        def __get__(self): return self.thisptr.njeval
        def __set__(self, int njeval): self.thisptr.njeval = njeval

    property stoich_reac:
        def __get__(self): return self.thisptr.stoich_reac
        def __set__(self, vector[vector[int]] stoich_reac): self.thisptr.stoich_reac = stoich_reac

    property stoich_prod:
        def __get__(self): return self.thisptr.stoich_prod
        def __set__(self, vector[vector[int]] stoich_prod): self.thisptr.stoich_prod = stoich_prod

    property k:
        def __get__(self): return self.thisptr.k
        def __set__(self, vector[double] k): self.thisptr.k = k

    property D:
        def __get__(self): return self.thisptr.D
        def __set__(self, vector[double] D): self.thisptr.D = D
