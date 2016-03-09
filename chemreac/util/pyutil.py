"""
chemreac.util.pyutil
--------------------

Utility functions used throughout chemreac.

"""
from __future__ import (absolute_import, division, print_function)


import sys
import numpy as np
import time


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


class progress(object):
    """ Print a progress bar of dots

    Parameters
    ----------
    iterable: iterable
        must have :attr:`__len__`
    output: fileobject
        default: sys.stdout
    proc_time: bool
        show process time (in seconds) passed at end of iteration.

    Examples
    --------
    >>> vals = list(range(7))
    >>> squares = []
    >>> for val in progress(vals, proc_time=False):
    ...     squares.append(val**2)
    ...
    7: .......
    >>> squares
    [0, 1, 4, 9, 16, 25, 36]

    """

    def __init__(self, iterable, output=None, proc_time=True):
        if proc_time is True:
            try:
                proc_time = time.process_time  # Py 3
            except AttributeError:
                proc_time = time.clock  # Py 2
        self._proc_time = proc_time
        self.iterable = iterable
        self.output = output
        self._cur_pos = 0
        self._len = len(iterable)
        if proc_time:
            self._t0 = self._proc_time()

    def __iter__(self):
        print("%d: " % self._len, file=self.output, end='')
        return self

    def next(self):
        if self._cur_pos >= self._len:
            print(' (%.3f s)\n' % (self._proc_time() - self._t0)
                  if self._proc_time else '\n', file=self.output, end='')
            raise StopIteration
        else:
            self._cur_pos += 1
            try:
                print('.', file=self.output, end='', flush=True)
            except TypeError:
                print('.', file=self.output, end='')
                if self.output is None or self.output is sys.stdout:
                    sys.stdout.flush()
            return self.iterable[self._cur_pos - 1]

    __next__ = next
