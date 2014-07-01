# -*- coding: utf-8 -*-

import numpy as np

def solver_linear_error(y, atol, rtol, logy=False):
    solver_err = np.abs(y*rtol) + atol
    if logy:
        res = np.exp(y - solver_err), np.exp(y + solver_err)
    else:
        res = y - solver_err, y + solver_err
    return np.array(res)
