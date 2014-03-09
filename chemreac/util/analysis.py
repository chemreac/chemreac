# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

ls=['-',':','--', '-.']
c='krgbycm'

def _plot_analysis(cb, labels, sys, tout, yout, indices, axes=None, titles=None,
                   limit=1e-10, logx=True):
    if axes == None:
        axes = [plt.subplot(len(indices),1,i+1) for i in range(len(indices))]
    else:
        assert len(axes) == len(indices)

    row_out = cb(sys, tout, yout, indices)
    for i, ax in enumerate(axes):
        ax.set_yscale('symlog', linthreshy=limit)
        if logx: ax.set_xscale('log')
        for j, lbl in enumerate(labels):
            if not np.any(np.abs(row_out[:,i,j]) > limit): continue
            ax.plot(tout, row_out[:,i,j], label=lbl, c=c[j%len(c)], ls=ls[j%len(ls)])
        ax.legend(loc='best', prop={'size': 11})
        if titles: ax.set_title(titles[i])
    return axes


def _get_jac_row(sys, tout, yout, indices):
    Jtmp = np.empty((sys.n*2+1, sys.n*sys.N), order='F')
    row_out = np.empty((yout.shape[0], len(indices), sys.n))
    for i, y in enumerate(yout):
        sys.banded_packed_jac_cmaj(tout[i], y, Jtmp)
        row_out[i,:,:] = Jtmp[indices,:]
    return row_out


def _get_per_func_out(sys, tout, yout, indices):
    out = np.empty((yout.shape[0], len(indices), sys.nr))
    for i, y in enumerate(yout):
        for j, si in enumerate(indices):
            sys.per_rxn_contrib_to_fi(tout[i], y, si, out[i, j, :])
    return out


def plot_jacobian(sys, tout, yout, substances, **kwargs):
    indices = [ri if isinstance(ri, int) else sys.names.index(ri) for ri in substances]
    print_names = sys.tex_names or sys.names
    _plot_analysis(_get_jac_row, print_names, sys, tout,
                   yout, indices, axes=axes,
                   titles=[print_names[i] for i in indices], **kwargs)


def plot_per_reaction_contribution(sys, tout, yout, substances, **kwargs):
    indices = [ri if isinstance(ri, int) else sys.names.index(ri)\
               for ri in substances]
    print_names = sys.tex_names or sys.names
    _plot_analysis(
        _get_per_func_out,
        ['R'+str(i) + ': ' + sys.to_Reaction(i).render(
            dict(zip(sys.names, print_names))) for \
         i in range(sys.nr)],
        sys, tout, yout, indices, axes=axes,
        titles=[print_names[i] for i in indices], **kwargs)
