#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import cPickle as pickle
except ImportError:
    import pickle

from collections import OrderedDict
import gzip
from itertools import product
import os

import numpy as np

from chemreac.util.pyutil import progress
from analytic_diffusion import integrate_rd

default_varied = OrderedDict([
    ('nstencil', [3, 5, 7]),  # categorical
    ('N', list(range(32, 385, 16))),  # continuous
    ('method', ['bdf', 'adams'])  # categorical
])

constant = dict(
    geom='f', nspecies=10,
    D=2e-3, t0=.25, tend=1.5, x0=0.0, xend=1.0, center=.5,
    nt=42, logt=False, logy=False, logx=False,
    random=False, p=0, a=2.7,
    linterpol=False, rinterpol=False, n_jac_diags=0, num_jacobian=False,
    atol=1e-8, rtol=1e-10,
    efield=True, random_seed=42, mobility=0.01,
    plot=False, savefig='None', verbose=False, yscale='linear',
    vline_limit=100, solver='sundials', iter_type='default',
    linear_solver='gmres', ilu_limit=1.0
)


def integrate(**kwargs):
    texecs = []
    nrepeat = 21
    ndrop = 7
    for i in range(nrepeat):
        tout, yout, info, rmsd_over_atol, sys, rmsd = integrate_rd(**kwargs)
        texecs.append(info['texec'])
    for _ in range(ndrop):
        texecs.pop(texecs.index(max(texecs)))
    texecs = np.array(texecs)
    info['texec'] = avg = np.average(texecs)
    info['dtexec'] = np.sqrt(np.sum((texecs - avg)**2/(len(texecs) - 1)))
    info['tout'] = tout
    # info['yout'] = yout
    info['rmsd'] = rmsd
    info['rmsd_over_atol'] = np.sqrt(np.sum(rmsd_over_atol**2))
    for k in default_varied:
        info[k] = kwargs[k]
    return info


def main(varied=None):
    if varied is None:
        varied = default_varied
    results = {
        'varied_keys': list(default_varied.keys()),
        'varied_values': list(default_varied.values())
    }
    all_params = list(product(*varied.values()))
    for params in progress(all_params):
        kwargs = constant.copy()
        kwargs.update(dict(zip(varied.keys(), params)))
        results[params] = integrate(**kwargs)
    basename = os.path.splitext(os.path.basename(__file__))[0]
    pickle.dump(results, gzip.open(basename+'.pkl', 'wb'))


if __name__ == '__main__':
    main()
