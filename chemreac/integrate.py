# -*- coding: utf-8 -*-

import time

import numpy as np

from scipy.integrate import ode

from chemreac import DENSE, BANDED, SPARSE

def run(sys, y0, t0, tend, nt, mode=None, log_time=False, **kwargs):
    assert y0.size == sys.n*sys.N

    defaults = {'name': 'vode', 'method': 'bdf', 'atol': 1e-12,
                'rtol': 1e-6, 'with_jacobian': True,
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

    if nt == 0: nt = 1024 # Good plotting density

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
        jout = np.zeros((sys.n*3+1, sys.n*sys.N), order='F')
    else:
        raise NotImplementedError

    def jac(t, y, *j_args):
        jac.neval += 1
        sys.banded_packed_jac_cmaj(t, y, jout)
        return jout
    jac.neval = 0

    runner = ode(f, jac=jac)
    runner.set_integrator(**kwargs)
    runner.set_initial_value(y0.flatten(), t0)
    if log_time:
        tout = np.logspace(t0+1e-12, np.log10(tend), nt+1)
    else:
        tout = np.linspace(t0, tend, nt+1)
    yout = np.empty((nt+1, sys.n*sys.N))
    yout[0,:] = y0
    texec = time.time()
    for i in range(1, nt+1):
        runner.integrate(tout[i])
        yout[i, :] = runner.y
    texec = time.time() - texec
    return tout, yout, {
        'texec': texec,
        'neval_f': f.neval,
        'neval_j': jac.neval,
    }
