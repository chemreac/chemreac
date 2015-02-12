# -*- coding: utf-8 -*-

"""
chemreac.util.grid
------------------
Grid related utilities for one dimensional grid of arbitrary spacing.
"""

from __future__ import print_function, division

import numpy as np


def padded_centers(x, nsidep):
    """
    Parameters
    ==========
    x: sequence
        strictly monotonically increasing sequence of positions of
        bin separators.
    nsidep: integer
        number of padding bins: (nstencil-1)/2
    """
    xc = x[:-1] + np.diff(x)/2
    return np.concatenate((
        2*x[0]-xc[:nsidep][::-1], xc, 2*x[-1]-xc[-nsidep:][::-1]
    ))


def pxci_to_bi(nstencil, N):
    """
    Generates a translation list converting x center indicesex starting
    at 0, which includes padding bins and into bin indices.

    Parameters
    ==========
    nstencil: integer
        Number of stencil points used
    N: integer
        Number of bins

    Returns
    =======
    list of bin indices.

    """
    nsidep = (nstencil-1)//2
    return list(range(nsidep)[::-1]) + list(range(N)) + list(
        range(N-1, N-nsidep-1, -1))


def stencil_pxci_lbounds(nstencil, N, lrefl=False, rrefl=False):
    """
    Generates a list of lower bounds in padded centers for each bin index
    for use in fintie difference scheme.

    Parameters
    ----------
    nstencil: int
        Number of stencil points used
    N: int
        Number of bins
    lrefl, rrefl: bool
        left and right reflective boundaries
    """
    nsidep = (nstencil-1)//2
    le = 0 if lrefl else nsidep
    re = 0 if rrefl else nsidep
    return [max(le, min(N + 2*nsidep - re - nstencil, i))
            for i in range(N)]
