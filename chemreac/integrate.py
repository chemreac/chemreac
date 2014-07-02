# -*- coding: utf-8 -*-

import time

import numpy as np

from scipy.integrate import ode

from chemreac import DENSE, BANDED, SPARSE


def run(sys, y0, tout, mode=None, **kwargs):
    """
    tout: at what times to report, e.g.:
        np.linspace(t0, tend, nt+1)
        np.logspace(t0+1e-12, np.log10(tend), nt+1)

    Returns
    =======
    yout: numpy array of shape (len(tout), sys.N, sys.n)
    """
    y0 = np.asarray(y0)
    assert y0.size == sys.n*sys.N

    defaults = {'name': 'vode', 'method': 'bdf', 'atol': 1e-9,
                'rtol': 1e-7, 'with_jacobian': True}

    if mode == None:
        if sys.N == 1:
            mode = DENSE
        elif sys.N > 1:
            mode = BANDED
        else:
            raise NotImplementedError

    if mode == BANDED:
        defaults['lband'] = sys.n
        defaults['uband'] = sys.n

    for k, v in defaults.items():
        if not k in kwargs:
            kwargs[k] = v

    # Create python callbacks with right signature
    fout = np.empty(sys.n*sys.N)
    def f(t, y, *f_args):
        # Python function closure circumvents reallocation
        f.neval += 1
        sys.f(t, y, fout)
        return fout
    f.neval = 0

    if mode == DENSE:
        jout = sys.alloc_jout(banded=False, order='F')
    elif mode == BANDED:
        # Currently SciPy <= v0.14 needs extra padding
        jout = sys.alloc_jout(banded=True, order='F', pad=sys.n)
    else:
        raise NotImplementedError

    def jac(t, y, *j_args):
        jac.neval += 1
        jout[...] = 0  # <--- this is very important (clear old LU decomp)
        if mode == DENSE:
            sys.dense_jac_cmaj(t, y, jout)
        else:
            sys.banded_packed_jac_cmaj(t, y, jout)
        return jout
    jac.neval = 0

    runner = ode(f, jac=jac if kwargs['with_jacobian'] else None)
    runner.set_integrator(**kwargs)
    runner.set_initial_value(y0.flatten(), tout[0])

    yout = np.empty((len(tout), sys.n*sys.N), order='C')
    yout[0,:] = y0
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
    return yout.reshape((len(tout), sys.N, sys.n)), info
