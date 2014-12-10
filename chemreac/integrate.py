# -*- coding: utf-8 -*-
"""
integrate
=========

This module provides functions for integrating the
system of ODEs which the ReactionDiffusion represent.
Currently the user may choose from using a
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)


import time

import numpy as np

from scipy import __version__ as __scipy_version__
from scipy.integrate import ode
scipy_version = tuple(map(int, __scipy_version__.split('.')[:3]))

from chemreac import DENSE, BANDED, SPARSE
from chemreac.util.analysis import suggest_t0


DEFAULTS = {
    'atol': 1e-9,
    'rtol': 1e-7,
}


def _integrate_cvode_direct(rd, y0, tout, mode=None, **kwargs):
    """
    see integrate.

    kwargs:
      method: linear multistep method: 'bdf' or 'adams'

    """
    from ._chemreac import cvode_direct

    # Handle kwargs
    new_kwargs = {}
    if mode is not None:
        raise NotImplementedError(
            "Sundials integrator auto-selects banded for N>1")
    atol = np.asarray(kwargs.pop('atol', DEFAULTS['atol']))
    if atol.ndim == 0:
        atol = atol.reshape((1,))
    new_kwargs['atol'] = atol
    new_kwargs['rtol'] = kwargs.pop('rtol', DEFAULTS['rtol'])
    new_kwargs['method'] = kwargs.pop('method', 'bdf')
    if kwargs.pop('with_jacobian', True) is False:
        raise ValueError("CVODE(S) wrapper uses explicit jacobian")
    if kwargs != {}:
        raise KeyError("Unkown kwargs: {}".format(kwargs))

    # Run the integration
    rd.neval_f = 0
    rd.neval_j = 0
    texec = time.time()
    try:
        yout = cvode_direct(rd, np.asarray(y0).flatten(),
                            np.asarray(tout).flatten(),
                            **new_kwargs)
    except RuntimeError:
        yout = np.ones((len(tout), rd.N, rd.n), order='C')/0  # NaN
        success = False
    else:
        success = True
    texec = time.time() - texec

    info = new_kwargs.copy()
    info.update({
        'neval_f': rd.neval_f,
        'neval_j': rd.neval_j,
        'texec': texec,
        'success': success
    })
    return yout, info


def _integrate_rk4(rd, y0, tout, **kwargs):
    """
    For demonstration purposes only, fixed step size
    give no error control and requires excessive work
    for accurate solutions. Unknown kwargs are simply
    ignored.

    see integrate
    """
    from ._chemreac import rk4
    texec = time.time()
    yout, Dyout = rk4(rd, y0, tout)
    texec = time.time() - texec
    info = {
        'neval_f': 4*(tout.size-1),
        'neval_j': 0,
        'texec': texec,
        'success': True,
    }
    return yout, info


def _integrate_scipy(rd, y0, tout, mode=None,
                     integrator_name='vode', **kwargs):
    """
    see integrate

    tout: array-like
        at what times to report, e.g.:
        - np.linspace(t0, tend, nt+1)
        - np.logspace(np.log10(t0 + 1e-12), np.log10(tend), nt+1)

    Returns
    =======
    yout: numpy array of shape (len(tout), rd.N, rd.n)

    """
    new_kwargs = {}
    y0 = np.asarray(y0)
    if y0.size != rd.n*rd.N:
        fmstr = "y0.size (={})not compatible with rd.n*rd.N (={})"
        raise ValueError(fmtstr.format(y0.size, rd.n*rd.N))

    if mode is None:
        if rd.N == 1:
            mode = DENSE
        elif rd.N > 1:
            mode = BANDED
        else:
            raise NotImplementedError

    if mode == BANDED:
        new_kwargs['lband'] = rd.n
        new_kwargs['uband'] = rd.n

    new_kwargs['atol'] = kwargs.pop('atol', DEFAULTS['atol'])
    new_kwargs['rtol'] = kwargs.pop('rtol', DEFAULTS['rtol'])
    new_kwargs['method'] = kwargs.pop('method', 'bdf')
    new_kwargs['with_jacobian'] = kwargs.pop('with_jacobian', True)
    if kwargs != {}:
        raise KeyError("Unkown kwargs: {}".format(kwargs))

    # Create python callbacks with right signature
    fout = np.empty(rd.n*rd.N)

    def f(t, y, *f_args):
        # Python function closure circumvents reallocation
        f.neval += 1
        rd.f(t, y, fout)
        return fout
    f.neval = 0

    if mode == DENSE:
        jout = rd.alloc_jout(banded=False, order='F')
    elif mode == BANDED:
        if scipy_version[0] <= 0 and scipy_version[1] <= 14:
            # Currently SciPy <= v0.14 needs extra padding
            jout = rd.alloc_jout(banded=True, order='F', pad=rd.n)
        else:
            # SciPy >= v0.15 need no extra padding
            jout = rd.alloc_jout(banded=True, order='F')
    else:
        raise NotImplementedError

    def jac(t, y, *j_args):
        jac.neval += 1
        jout[...] = 0  # <--- this is very important (clear old LU decomp)
        if mode == DENSE:
            rd.dense_jac_cmaj(t, y, jout)
        else:
            rd.banded_packed_jac_cmaj(t, y, jout)
        return jout
    jac.neval = 0

    runner = ode(f, jac=jac if new_kwargs['with_jacobian'] else None)
    runner.set_integrator(integrator_name, **new_kwargs)
    runner.set_initial_value(y0.flatten(), tout[0])

    yout = np.empty((len(tout), rd.n*rd.N), order='C')
    yout[0, :] = y0
    texec = time.time()
    for i in range(1, len(tout)):
        runner.integrate(tout[i])
        yout[i, :] = runner.y
    texec = time.time() - texec

    info = new_kwargs.copy()
    info.update({
        'integrator_name': integrator_name,
        'success': runner.successful(),
        'texec': texec,
        'neval_f': f.neval,
        'neval_j': jac.neval,
    })
    return yout.reshape((len(tout), rd.N, rd.n)), info


def sigm(x, lim=150., n=8):
    """
    Algebraic sigmoid to avoid overflow/underflow of 'double exp(double)'
    """
    return x/((x/lim)**n+1)**(1./n)


class Integration(object):
    """
    Model kinetcs by integrating system of ODEs using
    user specified solver.

    Parameters
    ----------
    solver: string
        "cvode_direct" or "scipy" where scipy uses VODE
        as the solver.
    rd: ReactionDiffusion instance
    C0: array
        initial concentrations (unscaled, untransformed)
    tout: array
        times for which to report solver results (untransformed)

    scaling: float
        Scale concentrations (and rate constants)
    sigm_damp: bool or tuple of (lim: float, n: int)
        conditionally damp C0 with an algebraic sigmoid when rd.logy == True.
        s(x) = x/((x/lim)**n+1)**(1./n)
        if sigm==True then `lim` and `n` are the default of `sigm()`
    C0_is_log: bool
        If True: passed values in C0 are taken to be the natural logarithm of
        initial concentrations. If False and rd.logy == True: a very small
        number is added to C0 to avoid applying log to zero (see `tiny`).
    tiny: float
        added to C0 when rd.logy==True and C0_is_log==False. Note that
        if you explicitly want to avoid adding tiny you need to set it
        to zero (e.g. when manually setting any C0==0 to some epsilon).
    (default: numpy.finfo(np.float64).tiny)

    **kwargs:
        mode: not supported by Sundials solver (current wrapper
              code auto selects banded for N>1 and uses dense
              mode for N==1)
        atol: float or sequence
            absolute tolerance of solution
        rtol: float
            relative tolerance of solution

    Attributes
    ----------
    Cout: array
        linear output concentrations
    yout: array
        output from solver: log(concentrations) if rd.logy == True
    info: dict
        Information from solver. Guaranteed to contain:
            - 'texec': execution time in seconds.
            - 'atol': float or array, absolute tolerance(s).
            - 'rtol': float, relative tolerance


    Methods
    -------
    _integrate()
        performs the integration, automatically called by __init__

    """

    _callbacks = {
        'cvode_direct': _integrate_cvode_direct,
        'scipy': _integrate_scipy,
        'rk4': _integrate_rk4,
    }

    def __init__(self, solver, rd, C0, tout, scaling=1.0,
                 sigm_damp=False, C0_is_log=False, tiny=None,
                 **kwargs):
        if solver not in self._callbacks:
            raise KeyError("Unknown solver %s" % solver)
        self.solver = solver
        self.rd = rd
        self.C0 = np.asarray(C0).flatten()
        self.tout = tout
        self.scaling = scaling
        self.sigm_damp = sigm_damp
        self.C0_is_log = C0_is_log
        self.tiny = tiny or np.finfo(np.float64).tiny
        self.kwargs = kwargs
        self.yout = None
        self.info = None
        self.Cout = None
        self._sanity_checks()
        self._integrate()

    def _sanity_checks(self):
        if not self.C0_is_log:
            if not np.all(self.C0 >= 0):
                raise ValueError("Negative concentrations encountered in C0")

    def _integrate(self):
        """
        Performs the integration by calling the callback chosen by
        self.solver. If rd.logy == True, a transformation of self.C0 to
        log(C0) will be performed before running the integration (the same
        is done for self.tout / rd.logt == True).

        After the integration is done the attributes `Cout`, `info` and `yout`
        are set. Cout is guaranteed to be linear concentrations (transformed
        from yout by calling exp if rd.logy==True) and yout is the unprocessed
        output from the solver.
        """
        # Possibly scale the concentrations
        if self.scaling != 1.0:
            C0 = self.scaling*self.C0
            ori_k = self.rd.k
            self.rd.k = [k*self.scaling**-(len(sa) - 1) for k, sa in
                         zip(self.rd.k, self.rd.stoich_actv)]
        else:
            C0 = self.C0

        # Transform initial concentrations
        if self.rd.logy:
            if not self.C0_is_log:
                C0 = np.log(C0 + self.tiny)

            if self.sigm_damp is True:
                y0 = sigm(C0)
            elif isinstance(self.sigm_damp, tuple):
                y0 = sigm(C0, *self.sigm_damp)
            else:
                y0 = C0
        else:
            if self.C0_is_log:
                if self.sigm_damp is True:
                    y0 = np.exp(sigm(C0))
                elif isinstance(self.sigm_damp, tuple):
                    y0 = np.exp(sigm(C0, *self.sigm_damp))
                else:
                    y0 = np.exp(C0)
            else:
                y0 = C0

        # Transform time
        if self.tout[0] == 0.0 and self.rd.logt:
            t0_set = True
            t0 = suggest_t0(self.rd, y0)
            t = np.log(self.tout + t0)  # conserve total time
        else:
            t0_set = False
            t = np.log(self.tout) if self.rd.logt else self.tout

        # Run the integration
        self.yout, self.info = self._callbacks[self.solver](
            self.rd, y0, t, **self.kwargs)
        self.internal_t = t
        self.info['t0_set'] = t0 if t0_set else False

        # Back-transform integration output into linear concentration
        self.Cout = np.exp(self.yout) if self.rd.logy else self.yout
        if self.scaling != 1.0:
            self.Cout /= self.scaling
            self.rd.k = ori_k


def run(*args, **kwargs):
    """
    ``run`` is provided for environment variable directed solver choice.

    Set ``CHEMREAC_SOLVER`` to indicate what integrator to
    use (default: "scipy").

    Set ``CHEMREAC_SOLVER_KWARGS`` to a string which can be eval'd to
    a python dictionary. e.g. "{'atol': 1e-4, 'rtol'=1e-7}"
    """
    import os
    environ_kwargs = os.getenv('CHEMREAC_SOLVER_KWARGS', None)
    if environ_kwargs:
        environ_kwargs = eval(environ_kwargs)
        if not isinstance(environ_kwargs, dict):
            fmtstr = "CHEMREAC_SOLVER_KWARGS not evaluated to a dictinary: {}"
            raise TypeError(fmtstr.format(environ_kwargs))
        kwargs.update(environ_kwargs)
    return Integration(
        os.getenv('CHEMREAC_SOLVER', 'scipy'), *args, **kwargs)
