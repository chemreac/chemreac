# -*- coding: utf-8 -*-
"""
chemreac.util.banded
--------------------

this module contains functions to deal with banded matrices.
"""
from __future__ import (absolute_import, division, print_function)


import numpy as np


def get_banded(A, n, N, n_jac_diags=1, order='C', padded=False):
    """ Turns a dense matrix into a banded one

    Turns a dense matrix (n·N) × (n·N) into a banded matrix
    including the diagonal and n super-diagonals and n sub-diagonals.

    Parameters
    ----------
    A: 2-dimensional square matrix
    n: int
        sub-block dimension
    N: int
        number of super-blocks
    n_jac_diags: int
        number of diagonals (default: 1)
    order: {'C', 'F'}, optional
        C- or Fortran-contiguous
    padded: bool, optional
        default: False, if True: A is padded with n rows along the top

    Raises
    ------
    ValueError on mismatch of A.shape and n*N
    """
    if A.shape != (n*N, n*N):
        raise ValueError("Shape of A != (n*N, n*N)")
    nouter = n * n_jac_diags
    B = np.zeros(((3 if padded else 2)*nouter + 1, n*N), order=order)
    for ri in range(n*N):
        for ci in range(max(0, ri-nouter), min(n*N, ri+nouter+1)):
            B[(2 if padded else 1)*nouter+ri-ci, ci] = A[ri, ci]
    return B


def get_jac_row_from_banded(J, rows, n, n_jac_diags=1):
    """
    Extracts rows from a banded matrix J

    Parameters
    ----------
    J: 2-dimensional array
        Source matrix with banded storage.
    rows: sequence
        indices of rows to extract
    n: integer
        row length
    n_jac_diags: integer
        number of diagonals (default: 1)
    """
    out = np.empty((len(rows), n))
    nouter = n * n_jac_diags
    for ri in rows:
        for ci in range(n):
            out[rows.index(ri), ci] = J[nouter+ri-ci, ci]
    return out


def get_dense(A, n, N, padded=False, n_jac_diags=1):
    out = np.zeros((n*N, n*N))
    nouter = n * n_jac_diags
    diag_offset = 2*nouter if padded else nouter
    for ri in range(n*N):
        for ci in range(max(0, ri-nouter), min(n*N, ri+nouter+1)):
            out[ri, ci] = A[diag_offset+ri-ci, ci]
    return out
