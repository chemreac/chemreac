# -*- coding: utf-8 -*-

import json


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
        'kerr': rd.kerr
    }
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
    kerr = data.pop('kerr', None)
    rd = RD(**data)
    if kerr: rd.kerr = kerr
    return rd
