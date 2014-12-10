
def dict_with_defaults(dct, *args):
    """
    Modifies a dictionary in-place by populating key/value pairs present
    in the default dictionaries which have no key in original dictionary `dct`.
    Returns same dictionary. The function is useful for passing along keyword
    argument dictionaries between functions.

    Parameters
    ----------
    dct: dict or None
    *args: dictionaries

    Returns
    -------
    dct: (possibly modified) input dictionary

    Examples
    --------
    >>> d = {1: None}
    >>> e = dict_with_defaults(d, {2: []})
    >>> d == {1: None, 2: []}
    True
    >>> d is e
    True
    >>> f = {'a': 1, 'b': 3}
    >>> f == dict_with_defaults({'a': 1}, {'b': 2, 'a': 7}, {'b': 3})
    True
    >>> g = {42: True, 'b': 3}
    >>> g == dict_with_defaults(None, {42: True, 'b': 2}, {'b': 3})
    True
    >>> dict_with_defaults(None)
    {}
    """
    try:
        ori_dct_keys = dct.keys()
    except AttributeError:
        dct = {}
        ori_dct_keys = []
    new_dct = {}
    for defaults in args:
        for k, v in defaults.items():
            if k not in ori_dct_keys:
                new_dct[k] = v
    dct.update(new_dct)
    return dct
