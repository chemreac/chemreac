
def set_dict_defaults_inplace(dct, *args):
    """
    Modifies a dictionary in-place by populating key/value pairs present in the
    default dictionaries which have no key in original dictionary `dct`. Useful
    for passing along keyword argument dictionaries between functions.

    Parameters
    ----------
    dct: dict
    *args: dictionaries

    Returns
    -------
    dct: (possibly modified) input dictionary

    Examples
    --------
    >>> d = {1: None}
    >>> set_dict_defaults_inplace(d, {2: []})
    >>> d == {1: None, 2: []}
    True
    >>> f = {'a': 1, 'b': 3}
    >>> g = {'a': 1}
    >>> set_dict_defaults_inplace(g, {'b': 2, 'a': 7}, {'b': 3})
    >>> f == g
    True
    >>> h = {42: True, 'b': 3}
    >>> i = {}
    >>> set_dict_defaults_inplace(i, {42: True, 'b': 2}, {'b': 3})
    >>> h == i
    True
    """
    ori_dct_keys = dct.keys()
    new_dct = {}
    for defaults in args:
        for k, v in defaults.items():
            if k not in ori_dct_keys:
                new_dct[k] = v
    dct.update(new_dct)
