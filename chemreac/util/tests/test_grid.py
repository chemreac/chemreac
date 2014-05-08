import numpy as np
from chemreac.util.grid import stencil_pxci_lbounds, padded_centers, pxci_to_bi

def test_stencil_pxci_lbounds():
    b = stencil_pxci_lbounds(5, 7)
    assert b == [2, 2, 2, 3, 4, 4, 4]

    b = stencil_pxci_lbounds(3, 5, lrefl=True)
    assert b == [0, 1, 2, 3, 3]

    b = stencil_pxci_lbounds(3, 5, rrefl=True)
    assert b == [1, 1, 2, 3, 4]

    b = stencil_pxci_lbounds(5, 7, lrefl=True, rrefl=True)
    assert b == [0, 1, 2, 3, 4, 5, 6]

def test_pxci_to_bi():
    yi = pxci_to_bi(5, 7)
    assert yi == [1, 0, 0, 1, 2, 3, 4, 5, 6, 6, 5]

def test_padded_centers():
    b = np.arange(10,21)
    xc_ = padded_centers(b, 2)
    print(xc_)
    assert np.allclose(xc_, np.linspace(8.5, 21.5, 14))
