import json

from chemreac import DENSE

def dump(rd, path):
    fh = open(path, 'wt')
    data = {
        'n': rd.n,
        'N': rd.N,
        'stoich_reac': rd.stoich_reac,
        'stoich_prod': rd.stoich_prod,
        'k': rd.k,
        'D': rd.D
    }
    json.dump(data, fh)

def load(path, RD=None, **kwargs):
    """
    Loads a ReactionDiffusion instance from implementation RD
    and overwrites with kwargs if passed
    """
    if not RD:
        from chemreac.cpp_chem_wrapper import \
            PyReactionDiffusion as ReactionDiffusion
        RD = ReactionDiffusion
    fh = open(path, 'rt')
    data = json.load(fh)
    data.update(kwargs)
    if not 'mode' in data: data['mode'] = DENSE
    return RD(**data)
