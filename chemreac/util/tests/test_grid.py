import numpy as np
from chemreac.util.grid import lbounds, padded_centers, y_indices

def test_lbounds():
    b = lbounds(5, 7)
    assert b == [2, 2, 2, 3, 4, 4, 4]

    b = lbounds(3, 5, lrefl=True)
    assert b == [0, 1, 2, 3, 3]

    b = lbounds(3, 5, rrefl=True)
    assert b == [1, 1, 2, 3, 4]

    b = lbounds(5, 7, lrefl=True, rrefl=True)
    assert b == [0, 1, 2, 3, 4, 5, 6]

def test_y_indices():
    yi = y_indices(5, 7)
    assert yi == [1, 0, 0, 1, 2, 3, 4, 5, 6, 6, 5]

def test_padded_centers():
    b = np.arange(10,21)
    xc_ = padded_centers(b, 2)
    print(xc_)
    assert np.allclose(xc_, np.linspace(8.5, 21.5, 14))
