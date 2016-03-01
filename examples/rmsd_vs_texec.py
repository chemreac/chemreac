#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import cPickle as pickle
except ImportError:
    import pickle

from itertools import product

import numpy as np

from analytic_diffusion import integrate_rd

varied = {
    'N': [32, 64, 128, 256],
    'nstencil': [3, 5, 7],
}

constant = dict(
    geom='f', nspecies=10,
    D=2e-3, t0=.25, tend=1.5, x0=0.0, xend=1.0, center=.5,
    nt=42, logt=False, logy=False, logx=False,
    random=False, p=0, a=2.7,
    linterpol=False, rinterpol=False, num_jacobian=False,
    method='bdf', atol=1e-8, rtol=1e-10,
    efield=True, random_seed=42, mobility=0.01,
    plot=False, savefig='None', verbose=False, yscale='linear',
    vline_limit=100,
)

def integrate(**kwargs):
    tout, yout, info, rmsd_over_atol, sys = integrate_rd(**kwargs)
    info['tout'] = tout
    info['yout'] = yout
    info['rmsd_over_atol'] = np.sqrt(np.sum(rmsd_over_atol**2))
    for k in varied:
        info[k] = kwargs[k]
    return info

    def main():
    results = {}
    for params in product(*varied.values()):
        kwargs = constant.copy()
        kwargs.update(dict(zip(varied.keys(), params)))
        results[params] = integrate(**kwargs)
    pickle.dump(results, open('analytic_diffusion_results.pkl', 'wb'))

if __name__ == '__main__':
    main()
