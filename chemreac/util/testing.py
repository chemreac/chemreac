# -*- coding: utf-8 -*-

"""
chemreac.util.testing
---------------------
Utility module for unit testing.
"""

from __future__ import (absolute_import, division, print_function)
import numpy as np
from chemreac.util.analysis import solver_linear_error


def check_rd_integration_run(cb, forgiveness=10, **kwargs):
    """
    Tests whether numerical solution is within error-bounds
    of reference solution. User provided callback "cb" needs
    to return a tuple of (yout, linCref, rd, info).
    """
    yout, linCref, rd, info = cb(**kwargs)
    assert info['success']
    for i in range(rd.n):
        try:
            atol = info['atol'][i]
        except:
            atol = info['atol']

        try:
            rtol = info['rtol'][i]
        except:
            rtol = info['rtol']
        lb, ub = solver_linear_error(yout[..., i], rtol, atol, rd.logy,
                                     scale_err=forgiveness, expb=rd.expb)
        assert np.all(lb < linCref[..., i])
        assert np.all(ub > linCref[..., i])


def spat_ave_rmsd_vs_time(Cout, Cref):
    """ Spatially averaged root mean square deviation versus time"""
    if not Cout.shape == Cref.shape:
        raise ValueError("Incompatible shapes: {}, {}".format(
            Cout.shape, Cref.shape))
    N = Cout.shape[1]
    err = Cout - Cref
    return np.sqrt(np.sum(err**2 / N, axis=1))
