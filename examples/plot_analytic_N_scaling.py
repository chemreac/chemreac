#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)


from chemreac import Geom_names
from collections import OrderedDict
import gzip
from itertools import chain, product
import os
import pickle

import numpy as np
import matplotlib.pyplot as plt

from chemreac.util.plotting import save_and_or_show_plot


def read(script_path, remove_prefix='plot_'):
    basename = os.path.splitext(os.path.basename(script_path))[0]
    if not basename.startswith(remove_prefix):
        raise ValueError("basename (%s) does not start with '%s'" %
                         (basename, remove_prefix))
    dirname = os.path.dirname(script_path)
    basename = basename[len('plot_'):]
    source = os.path.join(dirname, basename + '.pkl')
    if not os.path.exists(source):
        raise IOError("%s nonexistant. Run %s.py first" % (source, basename))
    results = pickle.load(gzip.open(source, 'rb'))
    varied_keys = results.pop('varied_keys')
    varied_vals = results.pop('varied_values')
    return results, OrderedDict(zip(varied_keys, varied_vals))


def to_array(results, varied, categorical, continuous, prop, cb=lambda x: x):
    def _len(key):
        return len(varied[key])
    arr = np.zeros(tuple(map(_len, categorical)) + (len(varied[continuous]),))

    def _index(pars):
        d = dict(zip(varied.keys(), pars))
        return tuple([varied[k].index(d[k]) for k
                      in chain(categorical, [continuous])])

    for params, data in results.items():
        arr[_index(params)] = cb(data[prop])
    return arr


def plot_with_matplotlib(
        savefig='None', dpi=300, errorbar=False, linx=False,
        categorical='geom,nspecies,nstencil', continuous='N',
        liny=False, ylims='None', nfits='7,6,4', prop='rmsd_over_atol',
        prop_label='RMSD/atol'):

    c = 'rbkg'
    m = 'osdx'

    categorical = categorical.split(',')
    results, varied = read(__file__)
    arr = to_array(results, varied, categorical, continuous,
                   prop, lambda x: np.average(x))

    if ylims != 'None':
        ylims = [[float(_) for _ in ylim.split(',')]
                 for ylim in ylims.split(';')]
    nfits = [int(_) for _ in nfits.split(',')]

    nrows, ncols, nseries = [len(varied[categorical[idx]]) for idx in range(3)]
    xdata = varied[continuous]
    plt.figure(figsize=(4*ncols, 10))
    for ri, ci in product(range(nrows), range(ncols)):
        geom_name = Geom_names[varied['geom'][ri]]
        ax = plt.subplot(nrows, ncols, ri*ncols + ci + 1)
        ax.set_xlabel(continuous)
        ax.set_ylabel(prop_label)
        ax.set_xscale('log', basex=2)
        ax.set_yscale('log', basey=2)
        if ci == 0:  # nspecies == 1
            plt.title('{}'.format(geom_name))
        else:
            plt.title('{}, {} decay(s)'.format(geom_name, ci))

        for si in range(nseries):
            ydata = arr[ri, ci, si]
            plt.plot(xdata, ydata, marker=m[si], ls='None', c=c[si])
            nfit = nfits[si]
            logx, logy = np.log(xdata[:nfit]), np.log(ydata[:nfit])
            pf = np.polyfit(logx, logy, 1)
            ax.plot(xdata[:nfit], np.exp(np.polyval(pf, logx)), ls='--',
                    c=c[si], label='%s: %s' % (str(varied[categorical[2]][si]),
                                               str(round(-pf[0], 1))))
        if ylims != 'None':
            ax.set_ylim(ylims[0])
        ax.legend(loc='upper right', prop={'size': 10}, numpoints=1)

    plt.tight_layout()
    save_and_or_show_plot(savefig=savefig, dpi=dpi)


if __name__ == '__main__':
    # e.g:
    #
    #  $ python plot_rmsd_vs_texec.py
    #  $ python plot_rmsd_vs_texec.py --savefig figure.png --dpi 100
    import argh
    argh.dispatch_command(plot_with_matplotlib)
