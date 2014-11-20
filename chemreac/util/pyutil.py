
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
    >>> d
    {1: None, 2: []}
    >>> d is e
    True
    >>> dict_with_defaults({'a': 1}, {'b': 2, 'a': 7}, {'b': 3})
    {'a': 1, 'b': 3}
    >>> dict_with_defaults(None, {42: True, 'b': 2}, {'b': 3})
    {42: True, 'b': 3}
    >>> dict_with_defaults(None)
    {}
    """
    dct = dct or {}
    ori_dct_keys = dct.keys()
    for defaults in args:
        for k, v in defaults.items():
            if k not in ori_dct_keys:
                dct[k] = v
    return dct
