# -*- coding: utf-8 -*-
"""
chemreac.integrate
==================

This module provides functions for integrating the
system of ODEs which the :py:class:`~chemreac.core.ReactionDiffusion` object
represents. The main class representing a numerical integration (for a set of
parameters) of the system of ODEs is :py:class:`Integration`.

If one does not want to hard code the choice of integrator and solver
parameters (e.g. tolerances), one may use :py:func:`run` which defers those
choices to the user of the script through the use of environment variables.

.. note :: Preferred ways to perform the integration is
    using :py:class:`Integration` or :py:func:`run`

"""

from __future__ import (absolute_import, division, print_function)


import os
import time

import numpy as np

from chemreac.units import get_derived_unit, to_unitless
from chemreac.util.analysis import suggest_t0


DEFAULTS = {
    'atol': 1e-9,
    'rtol': 1e-7,
}


class IntegrationError(Exception):
    pass


def integrate_cvode(rd, y0, tout, dense_output=None, **kwargs):
    """
    see :py:func:`integrate`

    kwargs:
      method: linear multistep method: 'bdf' or 'adams'

    """
    from ._chemreac import cvode_predefined, cvode_adaptive

    # Handle kwargs
    kwargs['atol'] = np.asarray(kwargs.pop('atol', DEFAULTS['atol']))
    if kwargs['atol'].ndim == 0:
        kwargs['atol'] = kwargs['atol'].reshape((1,))
    kwargs['rtol'] = kwargs.pop('rtol', DEFAULTS['rtol'])
    kwargs['method'] = kwargs.pop('method', 'bdf')
    kwargs['iter_type'] = {
        'default': 0, 'functional': 1, 'newton': 2}[
            kwargs.pop('iter_type', 'default').lower()]
    kwargs['linear_solver'] = {
        'default': 0, 'dense': 1, 'banded': 2, 'gmres': 10,
        'gmres_classic': 11, 'bicgstab': 20, 'tfqmr': 30}[
            kwargs.pop('linear_solver', 'default').lower()]
    if dense_output is None:
        dense_output = (len(tout) == 2)

    # Run the integration
    rd.zero_counters()
    time_wall = time.time()
    time_cpu = time.clock()
    try:
        if dense_output:
            if not len(tout) == 2:
                raise ValueError("dense_output implies tout == (t0, tend)")
            tout, yout = cvode_adaptive(
                rd, np.asarray(y0).flatten(), tout[0], tout[-1],
                kwargs.pop('atol'), kwargs.pop('rtol'), kwargs.pop('method'),
                **kwargs)
        else:
            yout = cvode_predefined(rd, np.asarray(y0).flatten(),
                                    np.asarray(tout).flatten(),
                                    **kwargs)
    except RuntimeError:
        yout = np.empty((len(tout), rd.N, rd.n), order='C')/0  # NaN
        success = False
    else:
        success = True
    time_wall = time.time() - time_wall
    time_cpu = time.clock() - time_cpu

    kwargs.update({
        'nfev': rd.nfev,
        'njev': rd.njev,
        'time_wall': time_wall,
        'time_cpu': time_cpu,
        'success': success
    })
    if kwargs['linear_solver'] >= 10:
        kwargs['nprec_setup'] = rd.nprec_setup
        kwargs['nprec_solve'] = rd.nprec_solve
        kwargs['njacvec_dot'] = rd.njacvec_dot
        kwargs['nprec_solve_ilu'] = rd.nprec_solve_ilu
        kwargs['nprec_solve_lu'] = rd.nprec_solve_lu
    kwargs.update(rd.last_integration_info)
    return yout, tout, kwargs


def _integrate_rk4(rd, y0, tout, **kwargs):
    """
    For demonstration purposes only, fixed step size
    give no error control and requires excessive work
    for accurate solutions. Unknown kwargs are simply
    ignored.

    see integrate
    """
    from ._chemreac import rk4
    time_wall = time.time()
    time_cpu = time.clock()
    yout, Dyout = rk4(rd, y0, tout)
    info = {
        'nfev': 4*(tout.size-1),
        'njev': 0,
        'time_wall': time.time() - time_wall,
        'time_cpu': time.clock() - time_cpu,
        'success': True,
    }
    return yout, tout, info


def _integrate_cb(callback, rd, y0, tout, linear_solver='dense',
                  dense_output=None, **kwargs):
    if dense_output is None:
        dense_output = (len(tout) == 2)
    if linear_solver != 'dense':
        raise NotImplementedError("Currently only dense jacobian is supported")
    new_kwargs = dict(y0=y0, dx0=1e-16*(tout[1]-tout[0]))
    new_kwargs.update(kwargs)
    if dense_output:
        new_kwargs['x0'] = tout[0]
        new_kwargs['xend'] = tout[1]
    else:
        new_kwargs['xout'] = tout
    info = {}
    info['atol'] = new_kwargs['atol'] = kwargs.pop('atol', DEFAULTS['atol'])
    info['rtol'] = new_kwargs['rtol'] = kwargs.pop('rtol', DEFAULTS['rtol'])

    def jac(t, y, jmat_out, dfdx_out):
        rd.dense_jac_rmaj(t, y, jmat_out)
        if rd.logt:
            fout = np.empty(rd.ny)
            rd.f(t, y, fout)
            dfdx_out[:] = fout
        else:
            dfdx_out[:] = 0
    new_kwargs['check_indexing'] = False
    time_wall = time.time()
    time_cpu = time.clock()
    if dense_output:
        xout, yout, info_ = callback[0](rd.f, jac, **new_kwargs)
    else:
        xout = tout
        yout, info_ = callback[1](rd.f, jac, **new_kwargs)
    info.update({
        'time_wall': time.time() - time_wall,
        'time_cpu': time.clock() - time_cpu,
        'success': True,
    })
    info.update(info_)
    return yout.reshape((xout.size, rd.N, rd.n)), xout, info


def _no_check(cb):
    def _cb(*args, **kwargs):
        kwargs['check_callable'] = False
        kwargs['check_indexing'] = False
        return cb(*args, **kwargs)
    return _cb


def integrate_pyodeint(*args, **kwargs):
    from pyodeint import integrate_adaptive, integrate_predefined
    return _integrate_cb((_no_check(integrate_adaptive),
                          _no_check(integrate_predefined)), *args, **kwargs)


def integrate_pygslodeiv2(*args, **kwargs):
    from pygslodeiv2 import integrate_adaptive, integrate_predefined
    return _integrate_cb((_no_check(integrate_adaptive),
                          _no_check(integrate_predefined)), *args, **kwargs)


def integrate_scipy(rd, y0, tout, linear_solver='default',
                    name='vode', dense_output=None,
                    **kwargs):
    """
    see :class:`Integration`

    Parameters
    ----------
    rd: ReactionDiffusion
    y0: array_like
        Initial conditions
    tout: array-like
        At what times to report, e.g.:
        - ``np.linspace(t0, tend, nt)``
        - ``np.logspace(np.log10(t0 + 1e-12), np.log10(tend), nt)``
    linear_solver: str (default: 'default')
        'dense' or 'banded'
    name: string (default: 'vode')
    dense_output: bool (default: None)
        if True, tout is taken to be length 2 tuple (t0, tend),
        if unspecified (None), length of tout decides (length 2 => True)

    Returns
    =======
    yout: numpy array of shape ``(len(tout), rd.N, rd.n)``.

    """

    from scipy import __version__ as __scipy_version__
    from scipy.integrate import ode
    scipy_version = tuple(map(int, __scipy_version__.split('.')[:2]))

    new_kwargs = {}
    y0 = np.asarray(y0)
    if y0.size != rd.n*rd.N:
        fmtstr = "y0.size (={})not compatible with rd.n*rd.N (={})"
        raise ValueError(fmtstr.format(y0.size, rd.n*rd.N))

    if linear_solver == 'default':
        if rd.N == 1:
            linear_solver = 'dense'
        elif rd.N > 1:
            linear_solver = 'banded'
    if linear_solver not in ('dense', 'banded'):
        raise NotImplementedError("Unkown linear_solver %s" % linear_solver)

    if linear_solver == 'banded':
        new_kwargs['lband'] = rd.n*rd.n_jac_diags
        new_kwargs['uband'] = rd.n*rd.n_jac_diags

    new_kwargs['atol'] = kwargs.pop('atol', DEFAULTS['atol'])
    new_kwargs['rtol'] = kwargs.pop('rtol', DEFAULTS['rtol'])
    new_kwargs['method'] = kwargs.pop('method', 'bdf')
    new_kwargs['with_jacobian'] = kwargs.pop('with_jacobian', True)
    new_kwargs['first_step'] = kwargs.pop('first_step', 0.0)
    if kwargs.pop('iter_type', 'default') != 'default':
        raise ValueError("iter_type unsupported by SciPy solver")
    if kwargs.pop('linear_solver', 'default') != 'default':
        raise ValueError("linear_solver unsupported by SciPy solver")
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

    if linear_solver == 'dense':
        jout = rd.alloc_jout(banded=False, order='F')
    elif linear_solver == 'banded':
        if scipy_version[0] <= 0 and scipy_version[1] <= 14:
            # Currently SciPy <= v0.14 needs extra padding
            jout = rd.alloc_jout(banded=True, order='F', pad=True)
        else:
            # SciPy >= v0.15 need no extra padding
            jout = rd.alloc_jout(banded=True, order='F')

    def jac(t, y, *j_args):
        jac.neval += 1
        jout[...] = 0  # <--- this is very important (clear old LU decomp)
        if linear_solver == 'dense':
            rd.dense_jac_cmaj(t, y, jout)
        else:
            if scipy_version[0] <= 0 and scipy_version[1] <= 14:
                # Currently SciPy <= v0.14 needs extra padding
                rd.banded_padded_jac_cmaj(t, y, jout)
            else:
                rd.banded_packed_jac_cmaj(t, y, jout)
        return jout
    jac.neval = 0

    runner = ode(f, jac=jac if new_kwargs['with_jacobian'] else None)
    runner.set_integrator(name, **new_kwargs)
    runner.set_initial_value(y0.flatten(), tout[0])

    if dense_output is None:
        dense_output = (len(tout) == 2)

    time_wall = time.time()
    time_cpu = time.clock()
    if dense_output:
        import warnings
        if not len(tout) == 2:
            raise ValueError("dense_output implies tout == (t0, tend)")
        # suppress warning printed by Fortran
        runner._integrator.iwork[2] = -1
        warnings.filterwarnings("ignore", category=UserWarning)
        yout = [y0]
        tstep = [tout[0]]
        while runner.t < tout[1]:
            runner.integrate(tout[1], step=True)
            tstep.append(runner.t)
            yout.append(runner.y)
        warnings.resetwarnings()
        tout = np.array(tstep)
        yout = np.array(yout)
    else:
        yout = np.empty((len(tout), rd.n*rd.N), order='C')
        yout[0, :] = y0
        for i in range(1, len(tout)):
            runner.integrate(tout[i])
            yout[i, :] = runner.y

    time_wall = time.time() - time_wall
    time_cpu = time.clock() - time_cpu

    info = new_kwargs.copy()
    info.update({
        'success': runner.successful(),
        'time_wall': time_wall,
        'time_cpu': time_cpu,
        'nfev': f.neval,
        'njev': jac.neval,
    })
    return yout.reshape((len(tout), rd.N, rd.n)), tout, info


def sigm(x, lim=150., n=8):
    r"""
    Algebraic sigmoid to avoid overflow/underflow of 'double exp(double)'.

    .. math ::

        s(x) = \frac{x}{\left((\frac{x}{lim})^n+1\right)^\frac{1}{n}}

    """
    return x/((x/lim)**n+1)**(1./n)


class Integration(object):
    """
    Model kinetcs by integrating system of ODEs using
    user specified integrator.

    Parameters
    ----------
    rd : ReactionDiffusion instance
    C0 : array
        Initial concentrations (untransformed, i.e. linear).
    tout : array
        Times for which to report results (untransformed).
    sigm_damp : bool or tuple of (lim: float, n: int)
        Conditionally damp C0 with an algebraic sigmoid when rd.logy == True.
        if sigm==True then `lim` and `n` are the default of :py:func:`sigm`.
    C0_is_log : bool
        If True: passed values in C0 are taken to be the natural logarithm of
        initial concentrations. If False and rd.logy == True: a very small
        number is added to C0 to avoid applying log to zero (see `tiny`).
    tiny : float
        Added to C0 when ``rd.logy==True`` and ``C0_is_log==False``. Note that
        if you explicitly want to avoid adding tiny you need to set it
        to zero (e.g. when manually setting any C0==0 to some epsilon).
        (default: None => ``numpy.finfo(np.float64).tiny``).
    integrator : string
        "cvode" or "scipy" where scipy uses VODE
        as the integrator.

    **kwargs :
        Keyword arguments passed on to integartor, e.g.:

        atol: float or sequence
            absolute tolerance of solution
        rtol: float
            relative tolerance of solution

    Attributes
    ----------
    Cout: array
        linear output concentrations
    yout: array
        output from integrator: log(concentrations) if rd.logy == True
    info: dict
        Information from integrator. Guaranteed to contain:
            - 'time_wall': execution time in seconds (wall clock).
            - 'time_cpu': execution time in seconds (cpu time).
            - 'atol': float or array, absolute tolerance(s).
            - 'rtol': float, relative tolerance
    rd: ReactionDiffusion instance
        same instance as passed in Parameters.

    Methods
    -------
    _integrate()
        performs the integration, automatically called by __init__


    """

    _callbacks = {
        'cvode': integrate_cvode,
        'scipy': integrate_scipy,
        'pyodeint': integrate_pyodeint,
        'pygslodeiv2': integrate_pygslodeiv2,
        'rk4': _integrate_rk4,
    }

    def __init__(self, rd, C0, tout, sigm_damp=False,
                 C0_is_log=False, tiny=None, integrator='scipy', **kwargs):
        if integrator not in self._callbacks:
            raise KeyError("Unknown integrator %s" % integrator)
        self.integrator = integrator
        self.rd = rd
        self.C0 = np.asarray(C0).flatten()
        self.tout = tout
        self.sigm_damp = sigm_damp
        self.C0_is_log = C0_is_log
        self.tiny = tiny or np.finfo(np.float64).tiny
        self.kwargs = kwargs
        self.yout = None
        self.info = None
        self.Cout = None
        self._sanity_checks()
        self._integrate()

    @classmethod
    def nondimensionalisation(cls, rd, C0, tout, **kw):
        if rd.unit_registry is None:
            raise ValueError("rd lacking unit_registry")

        def _n(arg, key):
            return to_unitless(arg, get_derived_unit(rd.unit_registry, key))
        return cls(rd, _n(C0, 'concentration'), _n(tout, 'time'), **kw)

    def _sanity_checks(self):
        if not self.C0_is_log:
            if np.any(self.C0 < 0):
                raise ValueError("Negative concentrations encountered in C0")

    def _integrate(self):
        """
        Performs the integration by calling the callback chosen by
        :attr:`integrator`. If rd.logy == True, a transformation of self.C0 to
        log(C0) will be performed before running the integration (the same
        is done for self.tout / rd.logt == True).

        After the integration is done the attributes `Cout`, `info` and `yout`
        are set. Cout is guaranteed to be linear concentrations (transformed
        from yout by calling exp if rd.logy==True) and yout is the unprocessed
        output from the integrator.
        """
        # Pre-processing
        # --------------
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
        tout = self.tout
        if tout[0] == 0.0 and self.rd.logt:
            t0_set = True
            t0 = suggest_t0(self.rd, y0)
            t = np.log(tout + t0)  # conserve total time
        else:
            t0_set = False
            t = np.log(tout) if self.rd.logt else tout

        # Run the integration
        # -------------------
        self.yout, self.internal_t, self.info = self._callbacks[
            self.integrator](self.rd, y0, t, **self.kwargs)
        self.info['t0_set'] = t0 if t0_set else False

        # Post processing
        # ---------------
        # Back-transform independent variable into linear time
        if self.rd.logt:
            self.tout = (np.exp(self.internal_t) - (t0 if t0_set else 0))
        else:
            self.tout = self.internal_t

        # Back-transform integration output into linear concentration
        self.Cout = np.exp(self.yout) if self.rd.logy else self.yout

    def get_with_units(self, attr):
        if attr == 'tout':
            return self.tout * get_derived_unit(self.rd.unit_registry, 'time')
        if attr == 'Cout':
            return self.Cout * get_derived_unit(self.rd.unit_registry,
                                                'concentration')
        else:
            raise ValueError("Unknown attr: %s" % attr)

    def internal_iter(self):
        """ Returns an iterator over (t, y)-pairs

        ``t`` is entries in ``internal_t`` and ``y`` is a (2-dim)
        array over the bins (1st dim) with the corresponding
        dependent variables (2nd dim)."""
        for idx, x in np.ndenumerate(self.internal_t):
            yield x, self.yout[idx, ...]


def run(*args, **kwargs):
    """
    ``run`` is provided for environment variable directed integration choices.

    Set ``CHEMREAC_INTEGRATION_KWARGS`` to a string which can be evaluated to
    a python dictionary. e.g. "{'integrator': 'cvode', 'atol': 1e-4}"
    """
    environ_kwargs = os.environ.get('CHEMREAC_INTEGRATION_KWARGS', None)
    if environ_kwargs:
        environ_kwargs = eval(environ_kwargs)
        if not isinstance(environ_kwargs, dict):
            fmtstr = "CHEMREAC_INTEGRATION_KWARGS not evaluated to a dict: {}"
            raise TypeError(fmtstr.format(environ_kwargs))
        kwargs.update(environ_kwargs)
    if kwargs.pop('nondimensionalisation', False):
        return Integration.nondimensionalisation(*args, **kwargs)
    else:
        return Integration(*args, **kwargs)
