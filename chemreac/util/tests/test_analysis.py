import numpy as np

from chemreac.util.analysis import solver_linear_error


def test_solver_linear_error():
    assert np.allclose(solver_linear_error(1.0, 1.0, 1.0), [-1, 3])
