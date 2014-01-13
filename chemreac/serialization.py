import json

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
        from . import ReactionDiffusion
        RD = ReactionDiffusion
    fh = open(path, 'rt')
    data = json.load(fh)
    data.update(kwargs)
    return RD(**data)
