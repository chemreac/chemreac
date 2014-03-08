# -*- coding: utf-8 -*-

import json

from . import ReactionDiffusion

def dump(rd, path):
    """
    Attributes ignored are:
    N, x, bin_k_factor, geom, logy, logt

    (geometrical factors and choice of variables)
    """
    fh = open(path, 'wt')
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
    }
    for attr in ReactionDiffusion.extra_attrs:
        data[attr] = getattr(rd, attr)
    json.dump(data, fh)


def load(path, RD=None, **kwargs):
    """
    Loads a ReactionDiffusion instance from implementation RD
    and overwrites with kwargs if passed
    """
    if not RD:
        from . import ReactionDiffusion
        RD = ReactionDiffusion
    fh = open(path, 'rt')
    data = json.load(fh)
    data.update(kwargs)
    extra_data = {}
    for attr in ReactionDiffusion.extra_attrs:
        extra_data[attr] = data.pop(attr, None)
    rd = RD(**data)
    for attr, val in extra_data.items():
        if val != None: setattr(rd, attr, val)
    return rd
