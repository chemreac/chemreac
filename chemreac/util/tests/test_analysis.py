import numpy as np
import pytest

from chemreac.util.analysis import (
    solver_linear_error, solver_linear_error_from_integration
)


def test_solver_linear_error():
    assert np.allclose(solver_linear_error(1.0, 1.0, 1.0), [-1, 3])


@pytest.mark.parametrize("atol", (1.0, [1.0]))
def test_solver_linear_error_from_integration(atol):
    class Dummy:
        unit_registry = None

        def expb(self, x):
            return np.exp(x)

    integr = Dummy()
    integr.yout = np.array([1.0]).reshape((1, 1, 1))
    print(integr.yout[slice(None), 0, 0])
    integr.info = {'rtol': 1.0, 'atol': atol}
    integr.rd = Dummy()
    integr.rd.logy = False
    assert np.allclose(solver_linear_error_from_integration(integr),
                       [[-1], [3]])
