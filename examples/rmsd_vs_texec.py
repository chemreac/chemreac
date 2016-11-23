#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


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
    vline_limit=100, integrator='cvode', iter_type='undecided',
    linear_solver='gmres', ilu_limit=1.0
)


def integrate(**kwargs):
    times_cpu, times_wal = [], []
    nrepeat = 21
    ndrop = 7
    for i in range(nrepeat):
        tout, yout, info, rmsd_over_atol, sys, rmsd = integrate_rd(**kwargs)
        times_cpu.append(info['time_cpu'])
        times_wal.append(info['time_wall'])
    for _ in range(ndrop):
        times_cpu.pop(times_cpu.index(max(times_cpu)))
        times_wal.pop(times_wal.index(max(times_wal)))
    times_cpu = np.array(times_cpu)
    times_wal = np.array(times_wal)
    info['time_cpu'] = avg_cpu = np.average(times_cpu)
    info['time_wall'] = avg_wal = np.average(times_wal)
    info['dtime_cpu'] = np.sqrt(
        np.sum((times_cpu - avg_cpu)**2 / (len(times_cpu) - 1)))
    info['dtime_wall'] = np.sqrt(
        np.sum((times_wal - avg_wal)**2 / (len(times_wal) - 1)))
    info['tout'] = tout
    # info['yout'] = yout
    info['rmsd'] = rmsd
    info['rmsd_over_atol'] = np.sqrt(np.sum(rmsd_over_atol**2))
    for k in default_varied:
        info[k] = kwargs[k]
    return info


def main(varied=None, verbose=False):
    if varied is None:
        varied = default_varied
    results = {
        'varied_keys': list(default_varied.keys()),
        'varied_values': list(default_varied.values())
    }
    all_params = list(product(*varied.values()))
    for params in progress(all_params) if verbose else all_params:
        kwargs = constant.copy()
        kwargs.update(dict(zip(varied.keys(), params)))
        results[params] = integrate(**kwargs)
    basename = os.path.splitext(os.path.basename(__file__))[0]
    pickle.dump(results, gzip.open(basename+'.pkl', 'wb'))


if __name__ == '__main__':
    try:
        import argh
        argh.dispatch_command(main)
    except ImportError:
        import sys
        if len(sys.argv) > 1:
            print("Unable to process parameters, argh missing. "
                  "Run 'pip install --user argh' to fix.", file=sys.stderr)
            sys.exit(os.EX_USAGE)  # non-ok exit
        main()
