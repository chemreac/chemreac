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

def load(path, RD=None, mode=0, N=None):
    if not RD:
        from chemreac.cpp_chem_wrapper import \
            PyReactionDiffusion as ReactionDiffusion
        RD = ReactionDiffusion
    fh = open(path, 'rt')
    data = json.load(fh)
    data['mode'] = mode
    if N: data['N'] = N
    return RD(**data)
