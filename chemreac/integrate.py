# -*- coding: utf-8 -*-

import time

import numpy as np

from scipy.integrate import ode

def run(sys, y0, t0, tend, nt, **kwargs):

    defaults = {'name': 'vode', 'method': 'bdf', 'atol': 1e-12,
                'rtol': 1e-6, 'with_jacobian': True,
                'first_step': 1e-9, 'lband': sys.n, 'uband': sys.n}

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

    jout = np.zeros((sys.n*3+1, sys.n*sys.N), order="F")
    def jac(t, y, *j_args):
        jac.neval += 1
        sys.banded_packed_jac_cmaj(t, y, jout)
        return jout
    jac.neval = 0

    # print('sys.stoich_reac', sys.stoich_reac)
    # print('sys.stoich_prod', sys.stoich_prod)
    # print('sys.stoich_actv', sys.stoich_actv)
    # print(y0)
    # f(t0, y0)
    # print(fout)

    runner = ode(f, jac=jac)
    runner.set_integrator(**kwargs)
    runner.set_initial_value(y0.flatten(), t0)
    tout = np.logspace(-9,np.log10(tend),nt)
    yout = np.empty((nt, sys.n*sys.N))
    texec = time.time()
    for i in range(nt):
        runner.integrate(tout[i])
        yout[i, :] = runner.y
    texec = time.time() - texec
    return tout, yout, {
        'texec': texec,
        'neval_f': f.neval,
        'neval_j': jac.neval,
    }
