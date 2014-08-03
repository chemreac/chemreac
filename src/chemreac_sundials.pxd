from libcpp.vector cimport vector
from chemreac cimport ReactionDiffusion

cdef extern from "chemreac_sundials.h" namespace "chemreac_sundials":
    cdef void direct_banded[T](ReactionDiffusion * sys,
                               const vector[T] atol,
                               const T rtol, const int lmm,
                               const T * const y0,
                               int nout,
                               const T * const tout,
                               T * const yout)
