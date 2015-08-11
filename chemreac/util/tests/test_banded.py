import numpy as np
import pytest

from chemreac.util.banded import get_banded, get_jac_row_from_banded, get_dense


def test_get_banded():
    A = np.array([[2.0, 3.0],
                  [5.0, 7.0]])
    B = get_banded(A, 1, 2)
    B_ref = np.array([
        [0.0, 3.0],
        [2.0, 7.0],
        [5.0, 0.0]
    ])
    assert np.allclose(B, B_ref)

    with pytest.raises(ValueError):
        get_banded(np.array([[1], [2]]), 1, 1)


def test_get_banded_F_padded():
    A = np.array([[2.0, 3.0],
                  [5.0, 7.0]])
    B = get_banded(A, 1, 2, order='F', padded=True)
    B_ref = np.array([
        [0.0, 0.0],
        [0.0, 3.0],
        [2.0, 7.0],
        [5.0, 0.0]
    ])
    assert np.allclose(B, B_ref)


def test_get_jac_row_from_banded():
    n = 3
    A = np.arange(n*n).reshape((n, n))
    B = get_banded(A, n, 1)
    assert np.allclose(get_jac_row_from_banded(B, [0], n), A[0, :])
    assert np.allclose(get_jac_row_from_banded(B, [1, 2], n), A[1:, :])


def test_get_dense():
    B_ref = np.array([
        [0.0, 3.0],
        [2.0, 7.0],
        [5.0, 0.0]
    ])
    assert np.allclose(get_dense(B_ref, 1, 2),  np.array([[2.0, 3.0],
                                                          [5.0, 7.0]]))
