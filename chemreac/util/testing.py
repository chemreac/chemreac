# -*- coding: utf-8 -*-

"""
testing
-------
Utility module for unit testing.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

import numpy as np


from chemreac.util.analysis import solver_linear_error


# _test... to avoid auto-detection by py.test
def _test_rd_integration_run(cb, forgiveness=10, **kwargs):
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
                                     scale_err=forgiveness)
        assert np.all(lb < linCref[..., i])
        assert np.all(ub > linCref[..., i])
