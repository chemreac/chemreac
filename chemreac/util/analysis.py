# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

ls=['-',':','--', '-.']
c='krgbycm'

def _plot_analysis(cb, labels, sys, tout, yout, indices, axes=None, titles=None,
                   limit=1e-10, logx=True, legend_kwargs=None):
    """
    pass legend_kwargs=False to supress legend
    """
    if legend_kwargs == None:
        legend_kwargs = dict(loc=None, prop={'size': 11})
    if axes == None:
        axes = [plt.subplot(len(indices),1,i+1) for i in range(len(indices))]
    else:
        assert len(axes) == len(indices)
    row_out = cb(sys, tout, yout, indices)
    for i, ax in enumerate(axes):
        ax.set_yscale('symlog', linthreshy=limit)
        if logx: ax.set_xscale('log')
        for j, lbl in enumerate(labels):
            if np.all(np.abs(row_out[:,i,j]) < limit): continue
            ax.plot(tout, row_out[:,i,j], label=lbl, c=c[j%len(c)], ls=ls[j%len(ls)])
        if legend_kwargs: ax.legend(**legend_kwargs)
        if titles: ax.set_title(titles[i])
    return axes


def _get_jac_row(sys, tout, yout, indices):
    Jtmp = np.zeros((sys.n*2+1, sys.n*sys.N), order='F')
    row_out = np.zeros((yout.shape[0], len(indices), sys.n))
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
                   yout, indices, titles=[print_names[i] for i in indices], **kwargs)


def plot_per_reaction_contribution(sys, tout, yout, substances, **kwargs):
    indices = [ri if isinstance(ri, int) else sys.names.index(ri)\
               for ri in substances]
    print_names = sys.tex_names or sys.names
    _plot_analysis(
        _get_per_func_out,
        [('*R' if i<np.sum(sys.bin_k_factor_span) else 'R')+str(i) + ': ' + \
         sys.to_Reaction(i).render(dict(zip(sys.names, print_names))) for \
         i in range(sys.nr)],
        sys, tout, yout, indices,
        titles=[print_names[i] for i in indices], **kwargs)


def _init_ax_substances_labels(sys, ax, substances, labels, xscale, yscale):
    # helper func..
    ax = ax or plt.subplot(1,1,1)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if substances == None:
        substances = range(sys.n)
    else:
        substances = [s if isinstance(s, int) else sys.names.index(s) for s in substances]

    if labels == None:
        labels = sys.tex_names or sys.names
        labels = [labels[i] for i in substances]
    else:
        assert len(labels) == len(substances)
    return ax, substances, labels


def plot_C_vs_t_in_bin(sys, tout, yout, bi, ax=None, labels=None, xscale='log',
                       yscale='log', substances=None,
                       titlefmt="C(t) in bin: {0:.3g} < x < {1:.3g}"):
    ax, substances, labels = _init_ax_substances_labels(sys, ax, substances,
                                                        labels, xscale, yscale)
    for i, lbl in zip(substances, labels):
        ax.plot(tout, yout[:, i+bi*sys.n], label=lbl, ls=ls[i%len(ls)], c=c[i%len(c)])
    ax.set_xlabel("t / s")
    ax.set_ylabel("C / M")
    ax.set_title(titlefmt.format(sys.x[bi], sys.x[bi+1]))
    ax.legend(loc='best', prop={'size': 11})
    return ax


def plot_C_vs_x(sys, tout, yout, substances, ti, ax=None, labels=None,
                xscale='log', yscale='log', basetitle="C(x)"):
    ax, substances, labels = _init_ax_substances_labels(sys, ax, substances,
                                                        labels, xscale, yscale)
    x_edges = np.repeat(sys.x, 2)[1:-1]
    for i, lbl in zip(substances, labels):
        y_edges = np.repeat(yout[ti, range(i, sys.n*sys.N, sys.n)], 2)
        ax.plot(x_edges, y_edges, label=lbl)
    ax.set_xlabel("x / m")
    ax.set_ylabel("C / M")
    ax.set_title(basetitle+" at t = {0:.3g} s".format(tout[ti]))
    ax.legend(loc='best', prop={'size': 11})


def plot_C_vs_t_and_x(sys, tout, yout, substance, ax=None, log10=False, **plot_kwargs):
    # it would be nice to accpet kwargs xscale='log', yscale='log', zscale='log'
    # but it's currently not supported by matplotlib:
    # http://matplotlib.1069221.n5.nabble.com/plot-surface-fails-with-log-axes-td10206.html
    substance = substance if isinstance(substance, int) else sys.names.index(substance)
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    ax = ax or plt.subplot(1,1,1, projection='3d')
    assert isinstance(ax, Axes3D)

    xty = [sys.x_centers, tout, yout]
    x_, t_, y_ = map(np.log10, xty) if log10 else xty
    X,T = np.meshgrid(x_, t_)
    if not 'cmap' in plot_kwargs: plot_kwargs['cmap'] = cm.gist_earth
    ax.plot_surface(X, T, y_[:, range(substance, sys.n*sys.N, sys.n)], **plot_kwargs)

    # fmtstr = "$\\mathrm{{\\log_{{10}}({})}}$" if log10 else "$\\mathrm{{{}}}$"
    fmtstr = "$log_{{10}}$({})" if log10 else "{}"
    ax.set_xlabel(fmtstr.format('x / m'))
    ax.set_ylabel(fmtstr.format('time / s'))
    ax.set_zlabel(fmtstr.format('C / M'))
    if sys.names:
        if sys.tex_names:
            name = sys.tex_names[substance]
        else:
            name = sys.names[substance]
        ax.set_title('['+name+'] vs. t and x')

    return ax

def plot_bin_k_factors(sys, ax=None, indices=None):
    ax = ax or plt.subplot(1,1,1)
    indices = indices or range(len(sys.bin_k_factor_span))
    factors = np.array(sys.bin_k_factor)
    x_edges = np.repeat(sys.x, 2)[1:-1]
    for i in indices:
        y_edges = np.repeat(factors[:,i], 2)
        ax.plot(x_edges, y_edges, label=i)
