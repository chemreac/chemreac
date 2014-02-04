# -*- coding: utf-8 -*-

__version__ = '0.0.4'

from ._chemreac import PyReactionDiffusion

DENSE, BANDED, SPARSE = range(3)
FLAT, SPHERICAL, CYLINDRICAL = range(3)

def ReactionDiffusion(
        n, stoich_reac, stoich_prod, k, N=0, D=None, x=None,
        stoich_actv=None, bin_k_factor=None, bin_k_factor_span=None,
        geom=FLAT, logy=False, logt=False):
    """
    Returns a PyReactionDiffusion instance (defined in _chemreac.pyx)

    Arguments:
    -`n`: number of species
    -`stoich_reac`: list of reactant index lists per reaction.
    -`stoich_prod`: list of product index lists per reaction.
    -`k`: array of reaction rate coefficients
    -`N`: number of compartments
    -`D`: diffusion coefficients (of length n)
    -`x`: compartment boundaries (of length N+1)
    -`stoich_actv`: list of ACTIVE reactant index lists per reaction.n
    -`bin_k_factor`: per compartment modulation of rate coefficients
    -`bin_k_factor_span`: spans over reactions affected by bin_k_factor
    -`geom`: any of (FLAT, SPHERICAL, CYLINDRICAL)

    The instance provides methods:

    f(t, y, fout)
    dense_jac_rmaj(t, y, Jout)
    dense_jac_cmaj(t, y, Jout)
    banded_jac_cmaj(t, y, Jout)
    banded_packed_jac_cmaj(t, y, Jout)

    some of which are used by chemreac.integrate.run
    """

    if N == 0:
        if x == None:
            N = 1
        else:
            N = len(x)-1

    if N > 1:
        assert n == len(D)
        _D = D
    else:
        _D = list([0]*n)

    if x == None: x = 1

    if isinstance(x, float) or isinstance(x, int):
        _x = [x/float(N)*i for i in range(N+1)]
    else:
        assert len(x) == N+1
        # monotonic:
        assert all([x[i+1]>x[i] for i in range(len(x)-1)])
        _x = x

    if stoich_actv == None:
        _stoich_actv = list([[]]*len(stoich_reac))
    else:
        _stoich_actv = stoich_actv
    assert len(_stoich_actv) == len(stoich_reac)

    assert len(stoich_reac) == len(stoich_prod) == len(k)
    assert geom in (FLAT, SPHERICAL, CYLINDRICAL)

    if geom == SPHERICAL or geom == CYLINDRICAL:
        assert _x[0] != 0.0

    # Handle bin_k_factor
    if bin_k_factor == None:
        assert bin_k_factor_span == None
        bin_k_factor = []
        bin_k_factor_span = []
    else:
        assert bin_k_factor_span != None
        assert len(bin_k_factor) == N
        assert all([len(x) == len(bin_k_factor_span) for x in bin_k_factor])
        assert all([x >= 0 for x in bin_k_factor_span])
    return PyReactionDiffusion(
        n, stoich_reac, stoich_prod, k, N, _D, _x, _stoich_actv,
        bin_k_factor, bin_k_factor_span, geom, logy, logt
    )
