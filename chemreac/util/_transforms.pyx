# -*- coding: utf-8 -*-

import cython
cimport numpy as cnp
import numpy as np

cdef extern double _sigmoid_algebraic_4 "sigmoid_algebraic_4" (double x)

@cython.wraparound(False)
@cython.boundscheck(False)
def sigmoid_algebraic_4_array(double[:] x):
    cdef int i
    cdef double v
    cdef cnp.ndarray[cnp.float64_t, ndim=1] out = np.empty(x.size, dtype=np.float64)
    for i in range(x.size):
        out[i] = _sigmoid_algebraic_4(x[i])
    return out

def sigmoid_algebraic_4(x):
    try:  # scalar
        return _sigmoid_algebraic_4(x)
    except:
        try:  # array like with flatten()
            return sigmoid_algebraic_4_array(x.flatten()).reshape(x.shape)
        except:  # x valid input to asarray
            arr = np.asarray(x, dtype=np.float64)
            return sigmoid_algebraic_4_array(arr.flatten()).reshape(arr.shape)
