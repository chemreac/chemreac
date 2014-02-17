# -*- coding: utf-8 -*-

import time

import numpy as np

from scipy.integrate import ode

from chemreac import DENSE, BANDED, SPARSE

def run(sys, y0, tout, mode=None, log_time=False, **kwargs):
    """
    tout: at what times to report, e.g.:
        np.linspace(t0, tend, nt+1)
        np.logspace(t0+1e-12, np.log10(tend), nt+1)
    """
    y0 = np.asarray(y0)
    assert y0.size == sys.n*sys.N

    defaults = {'name': 'vode', 'method': 'bdf', 'atol': 1e-12,
                'rtol': 1e-8, 'with_jacobian': True,
                'first_step': 1e-9}

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
        jout = np.zeros((sys.n*sys.N, sys.n*sys.N), order='F')
    elif mode == BANDED:
        jout = np.zeros((sys.n*3+1, sys.n*sys.N), order='F') # Currently SciPy need extra padding
    else:
        raise NotImplementedError

    def jac(t, y, *j_args):
        jac.neval += 1
        if mode == DENSE:
            sys.dense_jac_cmaj(t, y, jout)
        else:
            sys.banded_packed_jac_cmaj(t, y, jout)
        return jout
    jac.neval = 0

    runner = ode(f, jac=jac)
    runner.set_integrator(**kwargs)
    runner.set_initial_value(y0.flatten(), tout[0])
    yout = np.empty((len(tout), sys.n*sys.N))
    yout[0,:] = y0
    texec = time.time()
    for i in range(1, len(tout)):
        runner.integrate(tout[i])
        yout[i, :] = runner.y
    texec = time.time() - texec
    return yout, {
        'texec': texec,
        'neval_f': f.neval,
        'neval_j': jac.neval,
    }
