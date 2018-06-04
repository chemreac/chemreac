# -*- coding: utf-8 -*-
"""
chemreac.chemistry
==================

This module collects classes useful for describing substances,
reactions and reaction systems. The classes have methods to help
with consistent low-level conversion to numerical parameters of
the model. The classes are from the
`chempy <https://pypi.python.org/pypi/chempy>`_ package.

"""
from __future__ import print_function, division, absolute_import

from collections import OrderedDict
from chempy import Substance, Reaction, ReactionSystem


def mk_sn_dict_from_names(names, **kwargs):
    """
    Convenience function to generate a OrderedDict of Substance
    instances from a sequence of names and corresponding sequences
    of kwargs to Substance class.

    Parameters
    ----------
    names: sequence of strings
        names of substances
    **kwargs:
        sequences of corresponding keyword arguments

    Examples
    --------
    >>> d = mk_sn_dict_from_names(
    ...     'ABCD', D=[0.1, 0.2, 0.3, 0.4])
    >>> d  # doctest: +NORMALIZE_WHITESPACE
    OrderedDict([('A', <Substance(name=A, ...>),
    ('B', <Substance(name=B, ...)>), ('C', <Substance(name=C, ...)>),
    ('D', <Substance(name=D, ...)>)])
    >>> d['A'].name
    'A'
    """
    kwargs_list = []
    for i in range(len(names)):
        d = {}
        data = {}
        for k, v in kwargs.items():
            if k in Substance.attrs:
                d[k] = v[i]
            else:
                data[k] = v[i]
        d['data'] = data
        kwargs_list.append(d)

    return OrderedDict([(s, Substance(s, **kwargs_list[i])) for i, s
                        in enumerate(names)])
