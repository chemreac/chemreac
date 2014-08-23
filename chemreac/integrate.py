# -*- coding: utf-8 -*-

import time

import numpy as np

from scipy.integrate import ode

from chemreac import DENSE, BANDED, SPARSE


DEFAULTS = {
    'atol': 1e-9,
    'rtol': 1e-7,
}


def integrate_sundials(rd, y0, tout, mode=None, **kwargs):
    from ._chemreac import sundials_direct
    if mode is not None:
        raise NotImplementedError(
            "Sundials integrator auto-selectes banded for N>1")
    atol = np.asarray(kwargs.get('atol', DEFAULTS['atol']))
    if atol.ndim == 0:
        atol = atol.reshape((1,))
    rtol = np.asarray(kwargs.get('rtol', DEFAULTS['rtol']))
    lmm = kwargs.get('lmm', 'bdf')
    rd.neval_f = 0
    rd.neval_j = 0
    texec = time.time()
    try:
        yout = sundials_direct(rd, atol, rtol, lmm, y0, tout)
    except RuntimeError:
        yout = np.ones((len(tout), rd.n*rd.N), order='C')/0
        success = False
    else:
        success = True
    texec = time.time() - texec
    info = kwargs.copy()
    info.update({
        'neval_f': rd.neval_f,
        'neval_j': rd.neval_j,
        'texec': texec,
        'success': success
    })
    return yout, info


def integrate_scipy(rd, y0, tout, mode=None, **kwargs):
    """
    tout: at what times to report, e.g.:
        np.linspace(t0, tend, nt+1)
        np.logspace(t0+1e-12, np.log10(tend), nt+1)

    Returns
    =======
    yout: numpy array of shape (len(tout), rd.N, rd.n)
    """
    y0 = np.asarray(y0)
    assert y0.size == rd.n*rd.N

    defaults = DEFAULTS.copy()
    defaults.update({'name': 'vode', 'method': 'bdf', 'with_jacobian': True})

    if mode is not None:
        if rd.N == 1:
            mode = DENSE
        elif rd.N > 1:
            mode = BANDED
        else:
            raise NotImplementedError

    if mode == BANDED:
        defaults['lband'] = rd.n
        defaults['uband'] = rd.n

    for k, v in defaults.items():
        if k not in kwargs:
            kwargs[k] = v

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
        # Currently SciPy <= v0.14 needs extra padding
        jout = rd.alloc_jout(banded=True, order='F', pad=rd.n)
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

    runner = ode(f, jac=jac if kwargs['with_jacobian'] else None)
    runner.set_integrator(**kwargs)
    runner.set_initial_value(y0.flatten(), tout[0])

    yout = np.empty((len(tout), rd.n*rd.N), order='C')
    yout[0, :] = y0
    texec = time.time()
    for i in range(1, len(tout)):
        runner.integrate(tout[i])
        yout[i, :] = runner.y
    texec = time.time() - texec

    info = kwargs.copy()
    info.update({
        'success': runner.successful(),
        'texec': texec,
        'neval_f': f.neval,
        'neval_j': jac.neval,
    })
    return yout.reshape((len(tout), rd.N, rd.n)), info


# `run` is provided as backwards compability / environment chosen integrator
def run(*args, **kwargs):
    import os
    use_sundials = os.getenv('USE_SUNDIALS', '0')
    use_sundials = (use_sundials == '1') or (use_sundials.lower() == 'true')
    if use_sundials:
        return integrate_sundials(*args, **kwargs)
    else:
        return integrate_scipy(*args, **kwargs)
