# -*- coding: utf-8 -*-

import numpy as np
from chemreac.util.banded import get_jac_row_from_banded
import matplotlib.pyplot as plt


def spy(rd):
    # Spy
    from matplotlib import pyplot as plt
    b = rd.spy()
    plt.spy(b)
    plt.show()

def coloured_spy(A, cmap_name='gray', ax = None):
    from matplotlib import pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from matplotlib.cm import get_cmap

    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    plt.imshow(A, cmap=get_cmap(cmap_name), interpolation='none')
    ax = plt.gca()

    ya = ax.get_yaxis()
    ya.set_major_locator(MaxNLocator(integer=True))
    xa = ax.get_xaxis()
    xa.set_major_locator(MaxNLocator(integer=True))
    plt.colorbar()
    return ax


ls=['-',':','--', '-.']
c='krgbycm'


def _get_jac_row_over_t(rd, tout, yout, indices, bi=0):
    Jout = np.zeros((rd.n*2+1, rd.n*rd.N), order='F')
    row_out = np.zeros((yout.shape[0], len(indices), rd.n))
    for i, y in enumerate(yout):
        rd.banded_packed_jac_cmaj(tout[i], y.flatten(), Jout)
        Jtmp = Jout[:, bi*rd.n:(bi + 1)*rd.n]
        row_out[i,:,:] = get_jac_row_from_banded(Jtmp, indices, rd.n)
    return row_out


def _get_per_func_out(rd, tout, yout, indices):
    out = np.empty((yout.shape[0], len(indices), rd.nr))
    for i, y in enumerate(yout):
        flat_y = y.flatten()
        for j, si in enumerate(indices):
            rd.per_rxn_contrib_to_fi(tout[i], flat_y, si, out[i, j, :])
    return out


def _plot_analysis(cb, labels, rd, tout, yout, indices, axes=None,
                   titles=None, limit=1e-10, logx=True, legend_kwargs=None):
    """
    pass legend_kwargs=False to supress legend
    """
    if legend_kwargs == None:
        legend_kwargs = dict(loc=None, prop={'size': 11})
    if axes == None:
        axes = [plt.subplot(len(indices),1,i+1) for i in range(len(indices))]
    else:
        assert len(axes) == len(indices)
    row_out = cb(rd, tout, yout, indices)
    for i, ax in enumerate(axes):
        ax.set_yscale('symlog', linthreshy=limit)
        if logx: ax.set_xscale('log')
        for j, lbl in enumerate(labels):
            if np.all(np.abs(row_out[:,i,j]) < limit): continue
            ax.plot(tout, row_out[:,i,j], label=lbl, c=c[j%len(c)],
                    ls=ls[j%len(ls)])
        if legend_kwargs: ax.legend(**legend_kwargs)
        if titles: ax.set_title(titles[i])
        ax.set_xlabel("t / s")
    return axes


def plot_jacobian(rd, tout, yout, substances, **kwargs):
    indices = [ri if isinstance(ri, int) else rd.names.index(ri) for\
               ri in substances]
    print_names = rd.tex_names or rd.names
    axes = _plot_analysis(_get_jac_row_over_t, print_names, rd, tout,
                   yout, indices, titles=[print_names[i] for i in indices],
                          **kwargs)
    [ax.set_ylabel(
        "$\\frac{\\partial r_{tot}}{\partial C_i}~/~s^{-1}$")
     for ax in axes]



def plot_per_reaction_contribution(rd, tout, yout, substances, **kwargs):
    indices = [ri if isinstance(ri, int) else rd.names.index(ri)\
               for ri in substances]
    print_names = rd.tex_names or rd.names
    axes = _plot_analysis(
        _get_per_func_out,
        [('*R' if i<np.sum(rd.bin_k_factor_span) else 'R')+str(i) + ': ' + \
         rd.to_Reaction(i).render(dict(zip(rd.names, print_names))) for
         i in range(rd.nr)],
        rd, tout, yout, indices,
        titles=[print_names[i] for i in indices], **kwargs)
    [ax.set_ylabel("Reaction rate / $M\cdot s^{-1}$") for ax in axes]


def _init_ax_substances_labels(rd, ax, substances, labels, xscale, yscale):
    # helper func..
    ax = ax or plt.subplot(1,1,1)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if substances == None:
        substances = range(rd.n)
    else:
        substances = [s if isinstance(s, int) else rd.names.index(s) for \
                      s in substances]

    if labels == None:
        labels = rd.tex_names or rd.names
        labels = [labels[i] for i in substances]
    else:
        assert len(labels) == len(substances)
    return ax, substances, labels


def plot_C_vs_t_in_bin(
        rd, tout, yout, bi, ax=None, labels=None,
        xscale='log', yscale='log', substances=None,
        ttlfmt=(r"C(t) in bin: {0:.2g} m $\langle$"+\
                r" x $\langle$ {1:.2g} m"), legend_kwargs=None):
    if legend_kwargs == None:
        legend_kwargs = dict(loc='best', prop={'size': 11})

    ax, substances, labels = _init_ax_substances_labels(
        rd, ax, substances, labels, xscale, yscale)
    for i, lbl in zip(substances, labels):
        ax.plot(tout, yout[:, bi, i], label=lbl,
                ls=ls[i%len(ls)], c=c[i%len(c)])
    ax.set_xlabel("t / s")
    ax.set_ylabel("C / M")
    ax.set_title(ttlfmt.format(rd.x[bi], rd.x[bi+1]))
    ax.legend(**legend_kwargs)
    return ax


def plot_C_vs_x(rd, tout, yout, substances, ti, ax=None, labels=None,
                xscale='log', yscale='log', basetitle="C(x)"):
    ax, substances, labels = _init_ax_substances_labels(
        rd, ax, substances, labels, xscale, yscale)
    x_edges = np.repeat(rd.x, 2)[1:-1]
    for i, lbl in zip(substances, labels):
        y_edges = np.repeat(yout[ti, :, i], 2)
        ax.plot(x_edges, y_edges, label=lbl)
    ax.set_xlabel("x / m")
    ax.set_ylabel("C / M")
    ax.set_title(basetitle+" at t = {0:.3g} s".format(tout[ti]))
    ax.legend(loc='best', prop={'size': 11})


def plot_C_vs_t_and_x(rd, tout, yout, substance, ax=None, log10=False,
                      **plot_kwargs):
    # it would be nice to accpet kwargs xscale='log', yscale='log', zscale='log'
    # but it's currently not supported by matplotlib:
    # http://matplotlib.1069221.n5.nabble.com/plot-surface-fails-with-log-axes-td10206.html
    substance = substance if isinstance(substance, int) else rd.names.index(substance)
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    ax = ax or plt.subplot(1,1,1, projection='3d')
    assert isinstance(ax, Axes3D)

    xty = [rd.xcenters, tout, yout]
    x_, t_, y_ = map(np.log10, xty) if log10 else xty
    X,T = np.meshgrid(x_, t_)
    if not 'cmap' in plot_kwargs: plot_kwargs['cmap'] = cm.gist_earth
    ax.plot_surface(X, T, y_[:, range(substance, rd.n*rd.N, rd.n)], **plot_kwargs)

    # fmtstr = "$\\mathrm{{\\log_{{10}}({})}}$" if log10 else "$\\mathrm{{{}}}$"
    fmtstr = "$log_{{10}}$({})" if log10 else "{}"
    ax.set_xlabel(fmtstr.format('x / m'))
    ax.set_ylabel(fmtstr.format('time / s'))
    ax.set_zlabel(fmtstr.format('C / M'))
    if rd.names:
        if rd.tex_names:
            name = rd.tex_names[substance]
        else:
            name = rd.names[substance]
        ax.set_title('['+name+'] vs. t and x')

    return ax

def plot_bin_k_factors(rd, ax=None, indices=None, **kwargs):
    ax = ax or plt.axes(**kwargs)
    indices = indices or range(len(rd.bin_k_factor_span))
    factors = np.array(rd.bin_k_factor)
    x_edges = np.repeat(rd.x, 2)[1:]
    for i in indices:
        y_edges = np.pad(np.repeat(factors[:,i], 2), (0, 1), 'constant')
        ax.plot(x_edges, y_edges, label=i)
