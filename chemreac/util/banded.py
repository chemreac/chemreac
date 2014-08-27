import numpy as np

def get_banded(A, n, N):
    """
    Turns a dense matrix (n*N)*(n*N) into a banded one N
    blocks (n*n) + 1 sub and 1 super diagonal
    """
    B = np.zeros((2*n+1, n*N))
    for ri in range(n*N):
        for ci in range(max(0, ri-n), min(n*N, ri+n+1)):
            B[n+ri-ci, ci] = A[ri, ci]
    return B


def get_jac_row_from_banded(J, rows, n):
    """
    Extracts a rows from a banded matrix J

    Parameters
    ==========
    J: 2-dimensional array
        Source matrix with banded storage.
    rows: sequence
        indices of rows to extract
    n: integer
        row length
    """
    out = np.empty((len(rows),n))
    for ri in rows:
        for ci in range(n):
            out[rows.index(ri), ci] = J[n+ri-ci, ci]
    return out
