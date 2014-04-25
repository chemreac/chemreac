from chemreac.util.grid import bounds

def test_bounds():
    b = bounds(5, 7)
    assert b == [(0, 5), (0, 5), (0, 5), (1,6), (2,7), (2,7), (2,7)]
