import numpy as np
from chemreac.util.grid import (
    stencil_pxci_lbounds, padded_centers, pxci_to_bi, generate_grid
)


def test_stencil_pxci_lbounds():
    b = stencil_pxci_lbounds(5, 7)
    assert b == [2, 2, 2, 3, 4, 4, 4]

    b = stencil_pxci_lbounds(3, 5, lrefl=True)
    assert b == [0, 1, 2, 3, 3]

    b = stencil_pxci_lbounds(3, 5, rrefl=True)
    assert b == [1, 1, 2, 3, 4]

    b = stencil_pxci_lbounds(5, 7, lrefl=True, rrefl=True)
    assert b == [0, 1, 2, 3, 4, 5, 6]

    b = stencil_pxci_lbounds(3, 4, rrefl=True)
    # [y0, y0, y1, y2, y3, y3]
    assert b == [1, 1, 2, 3]


def test_pxci_to_bi():
    yi = pxci_to_bi(5, 7)
    assert yi == [1, 0, 0, 1, 2, 3, 4, 5, 6, 6, 5]

    yi = pxci_to_bi(3, 4)
    assert yi == [0, 0, 1, 2, 3, 3]


def test_padded_centers():
    b = [0, 2, 6]
    xc_ = padded_centers(b, 1)
    assert np.allclose(xc_, [-1, 1, 4, 8])

    b = np.arange(10, 21)
    xc_ = padded_centers(b, 2)
    print(xc_)
    assert np.allclose(xc_, np.linspace(8.5, 21.5, 14))


def test_generate_grid():
    g1 = generate_grid(0, 1, 2)
    assert np.allclose(g1, [0, 0.5, 1])

    g2 = generate_grid(0, 1, 2, random=True)
    assert g2[0] == 0
    assert g2[1] > 0 and g2[1] < 1
    assert g2[2] == 1
