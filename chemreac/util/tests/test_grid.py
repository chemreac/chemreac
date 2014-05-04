import numpy as np
from chemreac.util.grid import bounds, padded_centers, y_indices

def test_bounds():
    b = bounds(5, 7)
    assert b == [(2, 7), (2, 7), (2, 7), (3, 8), (4, 9), (4, 9), (4, 9)]

    b = bounds(3, 5, lrefl=True)
    assert b == [(0, 3), (1, 4), (2, 5), (3, 6), (3, 6)]

    b = bounds(3, 5, rrefl=True)
    assert b == [(1, 4), (1, 4), (2, 5), (3, 6), (4, 7)]

    b = bounds(5, 7, lrefl=True, rrefl=True)
    assert b == [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
                 (5, 10), (6, 11)]

def test_y_indices():
    yi = y_indices(5, 7)
    assert yi == [1, 0, 0, 1, 2, 3, 4, 5, 6, 6, 5]

def test_padded_centers():
    b = np.arange(10,21)
    xc_ = padded_centers(b, 2)
    print(xc_)
    assert np.allclose(xc_, np.linspace(8.5, 21.5, 14))
