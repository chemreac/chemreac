# -*- coding: utf-8 -*-

import json

import numpy as np

from . import ReactionDiffusion
from .core import get_unit
from .units import (
    unit_registry_to_human_readable,
    unit_registry_from_human_readable
)


def dump(rd, dest):
    """
    Serializes ReactionDiffusion instance to json format.

    Parameters
    ----------
    rd: ReactionDiffusion instance
    dest: file object or path string

    Notes
    -----
    Attributes ignored are:
    N, x, bin_k_factor, geom, logy, logt

    (geometrical factors and choice of variables)
    """

    if not hasattr(dest, 'write'):
        fh = open(dest, 'wt')
    else:
        fh = dest

    data = {
        'n': rd.n,
        'stoich_active': rd.stoich_active,
        'stoich_prod': rd.stoich_prod,
        'k': list(rd._k),
        'D': list(rd._D),
        'mobility': list(rd._mobility),
        # 'x': rd.x.tolist(),
        'stoich_inactv': rd.stoich_inactv,
        'unit_registry': unit_registry_to_human_readable(rd.unit_registry),
        'g_values': list(rd._g_values),
        'g_value_parents': rd.g_value_parents,
        # 'fields': fields
    }
    for attr in ReactionDiffusion.extra_attrs:
        data[attr] = getattr(rd, attr)
    json.dump(data, fh)


def load(source, RD=None, **kwargs):
    """
    Creates a ``RD`` instance (see below) from json serialized
    data.

    Parameters
    ----------
    source: file object or path string
    RD: subclass of ReactionDiffusion (default: ReactionDiffusion)
    \*\*kwargs
        override parameters in source with kwargs
    """
    RD = RD or ReactionDiffusion

    if not hasattr(source, 'read'):
        fh = open(source, 'rt')
    else:
        fh = source

    data = json.load(fh)
    unit_registry = data['unit_registry'] = unit_registry_from_human_readable(
        data.get('unit_registry', None))
    if 'D' in data:
        data['D'] = np.array(data['D'])*get_unit(unit_registry, 'diffusion')
    if 'mobility' in data:
        data['mobility'] = np.array(data['mobility'])*get_unit(
            unit_registry, 'electrical_mobility')
    if 'g_values' in data:
        data['g_values'] = [elem*g_unit for elem, g_unit in
                            zip(data['g_values'], RD.g_units(
                                unit_registry, data['g_value_parents']))]
    reaction_orders = map(len, data['stoich_active'])
    kunits = RD.k_units(unit_registry, reaction_orders)
    data['k'] = [kval*kunit for kval, kunit in zip(data['k'], kunits)]
    data.update(kwargs)
    extra_data = {}
    for attr in ReactionDiffusion.extra_attrs:
        extra_data[attr] = data.pop(attr, None)
    rd = RD(**data)
    for attr, val in extra_data.items():
        if val is not None:
            setattr(rd, attr, val)
    return rd
