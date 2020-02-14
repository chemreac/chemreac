# -*- coding: utf-8 -*-
"""
chemreac.util.analysis
----------------------

Functions to analyze numierc output from integration.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from ..units import get_derived_unit


def solver_linear_error(y, rtol, atol, logy=False, scale_err=1.0, expb=None):
    """
    Returns linear estimated error bounds from numerical integration

    Parameters
    ==========
    y : array_like
         Output from integration (before back-transformation is applied)
    rtol : float
         Relative tolerance
    atol : float
         Absolute tolerance
    logy : bool
         Is y from a run with logarithmic concentration?
    scale_err : float
         scale estimated error bounds (useful for testing)
    expb : callback
       exponential function in base b (e or 2 depending on ``use_log2``)

    Returns
    =======
    Array of shape (2, len(y)) with rows corresponding to lower and
    upper bounds around y.

    .. note:: Assumes maximum mangitude of error be: \
    :math:`\\boldsymbol{e}_{max} = \\|\\boldsymbol{y} \
    \\cdot \\mathrm{rtol}\\| + \\mathrm{atol}`
    """
    solver_err = scale_err*(np.abs(y*rtol) + atol)
    if logy:
        res = expb(y - solver_err), expb(y + solver_err)
    else:
        res = y - solver_err, y + solver_err
    return np.array(res)


def solver_linear_error_from_integration(integration, ti=slice(None), bi=0,
                                         si=0, **kwargs):
    """
    Convenience function wrapping :func:`solver_linear_error`

    Parameters
    ----------
    integration: Integration instance
    """
    try:
        atol_i = integration.info['atol'][si]
    except TypeError:
        atol_i = integration.info['atol']
    return solver_linear_error(
        integration.yout[ti, bi, si],
        integration.info['rtol'],
        atol_i,
        integration.rd.logy,
        expb=integration.rd.expb,
        **kwargs
    ) * get_derived_unit(integration.rd.unit_registry, 'concentration')


def suggest_t0(rd, y0, max_f=1.0):
    """
    Suggests an appropriate initial time,
    useful when logy==True and logt==True,
    If suggested t0 > 1, 1 is returned.

    Parameters
    ==========
    rd: ReactionDiffusion instance
         System at hand
    y0: sequence
         initial values
    max_f: float
         upper bound of absolute value for largest element in for the
         inital step.
    """
    fout = rd.alloc_fout()
    rd.f(0, np.asarray(y0), fout)
    fout_maxabs = np.max(np.abs(fout))
    if fout_maxabs < max_f:
        return 1.0
    else:
        return max_f/fout_maxabs


def eval_jacobian(rd, x, y):
    """ Evaluates the Jacobian matrix

    Parameters
    ----------
    rd: ReactionDiffusion instance
         System at hand
    x: float
         value of the independent variable
    y: array
         values of the dependent variables (you may need to transform)

    """
    banded = (rd.N > 1)
    jout = rd.alloc_jout(banded=banded, order='F')
    if banded:
        rd.banded_packed_jac_cmaj(x, y.flatten(), jout)
    else:
        rd.dense_jac_cmaj(x, y.flatten(), jout)
    return jout
