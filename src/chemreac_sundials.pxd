from libcpp.vector cimport vector
from chemreac cimport ReactionDiffusion
from libcpp cimport bool

cdef extern from "chemreac_sundials.h" namespace "chemreac_sundials":
    cdef void cvode_direct[T, U](
        U * sys,
        const vector[T] atol,
        const T rtol, const int lmm,
        const T * const y0,
        size_t nout,
        const T * const tout,
        T * const yout) except +
