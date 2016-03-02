# -*- coding: utf-8 -*-
# -*- mode: cython -*-
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "cvodes_wrapper.hpp" namespace "cvodes_wrapper":
    cdef void simple_integrate[T, U](
        U * sys,
        const vector[T] atol,
        const T rtol, const int lmm,
        const T * const y0,
        size_t nout,
        const T * const tout,
        T * const yout,
        bool with_jacobian,
        int, 
        int,
        int,
        double) except +
