"""
chemreac.util.pyutil
--------------------

Utility functions used throughout chemreac.

"""

import numpy as np


def monotonic(arr, positive=0, strict=False):
    """
    Check monotonicity of a serie

    Parameters
    ----------
    arr: array_like
        Array to be checked for monotonicity
    positive: -1, 0 or 1 (default: 0)
        -1: negative, 1: positive, 0: either
    strict: bool (default: False)
        Disallow zero difference between neighboring instances

    Examples
    --------
    >>> monotonic([0, 0, -1, -2])
    True
    >>> monotonic([0, 0, 1, 2], strict=True)
    False
    >>> monotonic([1, 2, 3], -1)
    False

    Returns
    -------
    bool
    """
    if positive not in (-1, 0, 1):
        raise ValueError("positive should be either -1, 0 or 1")
    delta = np.diff(arr)
    if positive in (0, 1):
        if strict:
            if np.all(delta > 0):
                return True
        else:
            if np.all(delta >= 0):
                return True
    if positive in (0, -1):
        if strict:
            if np.all(delta < 0):
                return True
        else:
            if np.all(delta <= 0):
                return True
    return False


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
