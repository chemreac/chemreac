from chemreac.util.pyutil import dict_with_defaults


def test_dict_with_defaults():
    assert dict_with_defaults(None, {}, {}) == {}
    assert dict_with_defaults(None, {1: 2}, {1: 4, 2: 8}) == {1: 4, 2: 8}
    d1 = {}
    d2 = dict_with_defaults(d1, {1: 2})
    assert d2 is d1
    assert d2 == d1
    assert dict_with_defaults(
        {1: 1, 2: 4}, {2: 8, 3: 27}, {3: 81, 4: 256}) == {
            1: 1, 2: 4, 3: 81, 4: 256}
