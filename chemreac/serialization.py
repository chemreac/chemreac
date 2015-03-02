# -*- coding: utf-8 -*-

import json

from . import ReactionDiffusion
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
        'stoich_reac': rd.stoich_reac,
        'stoich_prod': rd.stoich_prod,
        'k': rd.k.tolist(),
        'D': rd.D.tolist(),
        # 'x': rd.x.tolist(),
        'stoich_actv': rd.stoich_actv,
        # 'bin_k_factor': rd.bin_k_factor.,
        'bin_k_factor_span': rd.bin_k_factor_span.tolist(),
        'units': unit_registry_to_human_readable(rd.units)
    }
    for attr in ReactionDiffusion.extra_attrs:
        data[attr] = getattr(rd, attr)
    json.dump(data, fh)


def load(source, RD=None, **kwargs):
    """
    Creates a `RD` instance from json serialized
    data (where `RD` is an implementation of ReactionDiffusion)

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
    data.update(kwargs)
    extra_data = {}
    for attr in ReactionDiffusion.extra_attrs:
        extra_data[attr] = data.pop(attr, None)
    rd = RD(**data)
    for attr, val in extra_data.items():
        if val is not None:
            setattr(rd, attr, val)
    return rd
