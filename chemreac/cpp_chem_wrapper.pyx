# distutils: language = c++

from libcpp.vector cimport vector

cdef extern from "cpp_chem.hpp":
    cdef cppclass ReactionDiffusion:
        int n, N, nr, mode, geom, nfeval, njeval
        vector[vector[int]] stoich_reac
        vector[vector[int]] stoich_actv
        vector[vector[int]] stoich_prod
        vector[double] k
        vector[double] D
        vector[double] x

        ReactionDiffusion(int,
                          int,
                          vector[vector[int]],
                          vector[vector[int]],
                          vector[vector[int]],
                          vector[double],
                          vector[double],
                          vector[double],
                          int,
                          int) except +
        void f(double, double*, double*)
        void dense_jac_rmaj(double, double*, double*, int)
        void dense_jac_cmaj(double, double*, double*, int)
        void banded_jac_cmaj(double, double*, double*, int)
        void banded_packed_jac_cmaj(double, double*, double*, int)

DEF DENSE=0
DEF BANDED=1
DEF SPARSE=2

DEF FLAT=0
DEF SPHERICAL=1
DEF CYLINDRICAL=2

cdef class PyReactionDiffusion:
    cdef ReactionDiffusion *thisptr

    def __cinit__(self, int n, int N,
                  vector[vector[int]] stoich_reac,
                  vector[vector[int]] stoich_prod,
                  vector[double] k,
                  vector[double] D,
                  x = None,
                  stoich_actv = None,
                  int mode=0,
                  int geom=0,
              ):
        cdef vector[vector[int]] _stoich_actv
        cdef vector[double] _x

        if stoich_actv == None:
            _stoich_actv = list([[]]*len(stoich_reac))
        else:
            _stoich_actv = stoich_actv

        x = x or 1

        if isinstance(x, float) or isinstance(x, int):
            _x = [x/float(N)*i for i in range(N+1)]
        else:
            assert len(x) == N+1
            _x = x

        assert len(stoich_reac) == len(stoich_prod) == len(k)
        assert n == len(D)
        assert mode in range(3)
        print(_stoich_actv)
        print(_x)
        self.thisptr = new ReactionDiffusion(
            n, N, stoich_reac, stoich_prod, _stoich_actv, k,
            D, _x, mode, geom)

    def __dealloc__(self):
        del self.thisptr

    def f(self, double t, double [::1] y, double [::1] dydt):
        self.thisptr.f(t, &y[0], &dydt[0])

    def dense_jac_rmaj(self, double t, double [::1] y,
                       double [:, ::1] J):
        self.thisptr.dense_jac_rmaj(t, &y[0], &J[0,0], J.shape[1])

    def dense_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] J):
        self.thisptr.dense_jac_cmaj(t, &y[0], &J[0,0], J.shape[1])

    def banded_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] J):
        self.thisptr.banded_jac_cmaj(t, &y[0], &J[0,0], J.shape[0])

    def banded_packed_jac_cmaj(self, double t, double [::1] y,
                       double [::1, :] J):
        self.thisptr.banded_packed_jac_cmaj(
            t, &y[0], &J[0,0], J.shape[0])

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
