from __future__ import print_function, division

import numpy as np

def padded_centers(x, nsidep):
    xc = x[:-1] + np.diff(x)/2
    return np.concatenate((
        2*x[0]-xc[:nsidep][::-1], xc, 2*x[-1]-xc[-nsidep:][::-1]
    ))

def y_indices(nstencil, N):
    """ Corresponding y indices for padded centers for refl """
    nsidep = (nstencil-1)//2
    return range(nsidep)[::-1]+range(N)+range(N-1, N-nsidep-1, -1)

def lbounds(nstencil, N, lrefl=False, rrefl=False):
    """
    Lower bounds in (possibly padded) centers for use in fintie difference.
    """
    nsidep = (nstencil-1)//2
    le = 0 if lrefl else nsidep
    re = 0 if rrefl else nsidep
    return [max(le, min(N + 2*nsidep - re - nstencil, i)) for i in range(N)]
