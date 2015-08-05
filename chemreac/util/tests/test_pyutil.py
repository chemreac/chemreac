from chemreac.util.pyutil import set_dict_defaults_inplace, monotonic


def test_set_dict_defaults_inplace():
    d1 = {}
    set_dict_defaults_inplace(d1, {}, {})
    assert d1 == {}

    d2 = {}
    set_dict_defaults_inplace(d2, {1: 2}, {1: 4, 2: 8})
    assert d2 == {1: 4, 2: 8}

    d3 = {}
    set_dict_defaults_inplace(d3, {1: 2})
    assert d3 == {1: 2}

    d4 = {1: 1, 2: 4}
    set_dict_defaults_inplace(
        d4, {2: 8, 3: 27}, {3: 81, 4: 256})
    assert d4 == {1: 1, 2: 4, 3: 81, 4: 256}


def test_monotonic():
    assert monotonic([0, 1, 2])
    assert monotonic([0, -1, -2])
    assert not monotonic([0, -1, -2, 3])
    assert monotonic([0, -1, -2, -3], strict=True)
    assert not monotonic([0, 1, 2, 3, 3], strict=True)
    assert not monotonic([0, -1, -2], 1)
    assert monotonic([0, -1, -2, -3], -1)
