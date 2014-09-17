# -*- coding: utf-8 -*-
"""
banded
------

Functions to deal with banded matrices.
"""

import numpy as np


def get_banded(A, n, N):
    """
    Turns a dense matrix (n*N)*(n*N) into a banded matrix
    including the diagonal and n super-diagonals and n sub-diagonals

    Parameters
    ----------
    A: 2-dimensional square matrix
    n: int
        sub-block dimension
    N: int
        number of super-blocks

    Raises
    ------
    ValueError on mismatch of A.shape and n*N
    """
    if A.shape != (n*N, n*N):
        raise ValueError("Shape of A != (n*N, n*N)")
    B = np.zeros((2*n+1, n*N))
    for ri in range(n*N):
        for ci in range(max(0, ri-n), min(n*N, ri+n+1)):
            B[n+ri-ci, ci] = A[ri, ci]
    return B


def get_jac_row_from_banded(J, rows, n):
    """
    Extracts a rows from a banded matrix J

    Parameters
    ----------
    J: 2-dimensional array
        Source matrix with banded storage.
    rows: sequence
        indices of rows to extract
    n: integer
        row length
    """
    out = np.empty((len(rows), n))
    for ri in rows:
        for ci in range(n):
            out[rows.index(ri), ci] = J[n+ri-ci, ci]
    return out
