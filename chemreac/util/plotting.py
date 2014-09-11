# -*- coding: utf-8 -*-

"""
plotting
--------

convenience functions to create matplotlib plots
of results.
"""

import numpy as np
from chemreac.util.banded import get_jac_row_from_banded
import matplotlib.pyplot as plt


def coloured_spy(A, cmap_name='gray', ax=None):
    """
    Convenience function for using matplotlib to
    generate a spy plot for inspecting e.g. jacobian or
    its LU decomposition.

    Parameters
    ----------
    A: 2D array
        Array to inspect, populated e.g. by jacobian callback.
    cmap_name: string (default: gray)
        name of matplotlib colormap to use
    ax: Axes instance (default: None)
         Axes to plot to, if not given a new figure will be generated.

    Returns
    -------
    Axes instance plotted to
    """
    from matplotlib.ticker import MaxNLocator
    from matplotlib.cm import get_cmap

    if ax is None:
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


DEFAULT = dict(
    ls=['-', ':', '--', '-.'],
    c='krgbycm'
)


def _get_jac_row_over_t(rd, tout, yout, indices, bi=0):
    Jout = np.zeros((rd.n*2+1, rd.n*rd.N), order='F')
    row_out = np.zeros((yout.shape[0], len(indices), rd.n))
    for i, y in enumerate(yout):
        rd.banded_packed_jac_cmaj(tout[i], y.flatten(), Jout)
        Jtmp = Jout[:, bi*rd.n:(bi + 1)*rd.n]
        row_out[i, :, :] = get_jac_row_from_banded(Jtmp, indices, rd.n)
    return row_out


def _get_per_func_out(rd, tout, yout, indices):
    out = np.empty((yout.shape[0], len(indices), rd.nr))
    for i, y in enumerate(yout):
        flat_y = y.flatten()
        for j, si in enumerate(indices):
            rd.per_rxn_contrib_to_fi(tout[i], flat_y, si, out[i, j, :])
    return out


def _plot_analysis(cb, labels, rd, tout, yout, indices, axes=None,
                   titles=None, lintreshy=1e-10, logx=True,
                   legend_kwargs=None, ls=None, c=None):
    """
    Parameters
    ----------
    cb: callback
        callback with signature (rd, tout, yout, indices) returning
        3-dimensional array with shape (tout.size, len(axes), len(labels))
    labels: sequence of strings
        labels of individual plots
    rd: ReactionDiffusion instance
    tout: 1-dimensional array of floats
    yout: solver output
    indices: 4 arguments for callback
    axes: sequence of matplotlib Axes instances
        (default: len(indices) number of subplot axes)
    titles: titles per axes
    lintreshy: float
        symlog option 'linthreshy' (default: 1e-10)
    logx: set x scale to 'log'
    legend_kwargs: dict
        dict of kwargs to pass to matplotlib legend function
        (default: {'loc': None, 'prop': {'size': 11}}), set
        to False to suppress legend.
    ls: sequence of strings
        linestyles
    c: sequence of strings
        colors
    """
    if legend_kwargs is None:
        legend_kwargs = dict(loc=None, prop={'size': 11})
    if axes is None:
        axes = [plt.subplot(len(indices), 1, i+1) for i in range(
            len(indices))]
    else:
        assert len(axes) == len(indices)
    ls = ls or DEFAULT['ls']
    c = c or DEFAULT['c']
    row_out = cb(rd, tout, yout, indices)
    for i, ax in enumerate(axes):
        ax.set_yscale('symlog', linthreshy=lintreshy)
        if logx:
            ax.set_xscale('log')
        for j, lbl in enumerate(labels):
            if np.all(np.abs(row_out[:, i, j]) < lintreshy):
                continue
            ax.plot(tout, row_out[:, i, j], label=lbl, c=c[j % len(c)],
                    ls=ls[j % len(ls)])
        if legend_kwargs is not False:
            ax.legend(**legend_kwargs)
        if titles:
            ax.set_title(titles[i])
        ax.set_xlabel("t / s")
    return axes


def plot_jacobian(rd, tout, yout, substances, **kwargs):
    """
    Plots time evolution of Jacobian values (useful when investigating
    numerical instabilities).

    Parameters
    ----------
    rd: ReactionDiffusion
        system at hand
    tout: array_like
        output time from integration
    yout: array_like
        output data from integration
    substances: iterable of int or string
        indices or names of substances to plot jacobian values for
    """
    indices = [si if isinstance(si, int) else rd.substance_names.index(si) for
               si in substances]
    print_names = rd.tex_names or rd.substance_names
    axes = _plot_analysis(_get_jac_row_over_t, print_names, rd, tout, yout,
                          indices, titles=[
                              print_names[i] for i in indices], **kwargs)
    [ax.set_ylabel(
        "$\\frac{\\partial r_{tot}}{\partial C_i}~/~s^{-1}$")
     for ax in axes]


def plot_per_reaction_contribution(rd, tout, yout, substances, **kwargs):
    """
    Plots contributions to concentration derivatives of selected
    substances from individual reactions.

    Parameters
    ----------
    rd: ReactionDiffusion
    tout: 1D array of floats
    yout: output from solver
    substances: sequence of Substance instances
    **kwargs: kwargs passed on to _plot_analysis

    Returns
    -------
    list of matplotlib.axes.Axes instances
    """
    indices = [ri if isinstance(ri, int) else rd.substance_names.index(ri)
               for ri in substances]
    print_names = rd.tex_names or rd.substance_names
    axes = _plot_analysis(
        _get_per_func_out,
        [('*R' if i < np.sum(rd.bin_k_factor_span) else 'R') +
         str(i) + ': ' + rd.to_Reaction(i).render(dict(zip(
             rd.substance_names, print_names))) for i in range(rd.nr)],
        rd, tout, yout, indices,
        titles=[print_names[i] for i in indices], **kwargs)
    for ax in axes:
        ax.set_ylabel(r"Reaction rate / $M\cdot s^{-1}$")
    return axes


def _init_ax_substances_labels(rd, ax, substances, labels, xscale, yscale):
    # helper func..
    ax = ax or plt.subplot(1, 1, 1)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if substances is None:
        substances = range(rd.n)
    else:
        substances = [s if isinstance(s, int) else rd.substance_names.index(s)
                      for s in substances]

    if labels is None:
        labels = rd.tex_names or rd.substance_names
        labels = [labels[i] for i in substances]
    else:
        assert len(labels) == len(substances)
    return ax, substances, labels


def plot_C_vs_t_in_bin(
        rd, tout, yout, bi, ax=None, labels=None,
        xscale='log', yscale='log', substances=None,
        ttlfmt=(r"C(t) in bin: {0:.2g} m $\langle$" +
                r" x $\langle$ {1:.2g} m"), legend_kwargs=None,
        ls=None, c=None):
    """
    Plots bin local concentration as function of time for selected
    substances.

    Parameters
    ----------
    rd: ReactionDiffusion
    tout: 1D array of floats
    yout: output from solver
    bi: bin index
    ax: Axes instance
    labels: sequence of strings
    xscale: matplotlib scale choice (e.g. 'log', 'symlog')
    yscale: matplotlib scale choice (e.g. 'log', 'symlog')
    substances: sequence of indies or names of substances
    ttlfmt: string formatted with bin boundaries (set to empty to suppress)
    legend_kwargs: dict
        kwargs passed to matplotlib legend function,
        (default: {'loc': None, 'prop': {'size': 11}}), set
        to False to suppress legend.
    ls: sequence of strings
        linestyles
    c: sequence of strings
        colors

    Returns
    =======
    Axes instance
    """
    if legend_kwargs is None:
        legend_kwargs = dict(loc='best', prop={'size': 11})
    ls = ls or DEFAULT['ls']
    c = c or DEFAULT['c']
    ax, substances, labels = _init_ax_substances_labels(
        rd, ax, substances, labels, xscale, yscale)
    for i, lbl in zip(substances, labels):
        ax.plot(tout, yout[:, bi, i], label=lbl,
                ls=ls[i % len(ls)], c=c[i % len(c)])
    ax.set_xlabel("t / s")
    ax.set_ylabel("C / M")
    if ttlfmt:
        ax.set_title(ttlfmt.format(rd.x[bi], rd.x[bi+1]))
    if legend_kwargs is not False:
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
    # it would be nice to accpet kwargs
    #    xscale='log', yscale='log', zscale='log'
    # but it's currently not supported by matplotlib:
    # http://matplotlib.1069221.n5.nabble.com/plot-surface-fails-with-log-axes-td10206.html
    substance = (substance if isinstance(substance, int) else
                 rd.substance_names.index(substance))
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    ax = ax or plt.subplot(1, 1, 1, projection='3d')
    assert isinstance(ax, Axes3D)

    xty = [rd.xcenters, tout, yout]
    x_, t_, y_ = map(np.log10, xty) if log10 else xty
    X, T = np.meshgrid(x_, t_)
    if 'cmap' not in plot_kwargs:
        plot_kwargs['cmap'] = cm.gist_earth
    ax.plot_surface(X, T, y_[:, range(substance, rd.n*rd.N, rd.n)],
                    **plot_kwargs)

    fmtstr = "$log_{{10}}$({})" if log10 else "{}"
    ax.set_xlabel(fmtstr.format('x / m'))
    ax.set_ylabel(fmtstr.format('time / s'))
    ax.set_zlabel(fmtstr.format('C / M'))
    if rd.substance_names:
        if rd.tex_names:
            name = rd.tex_names[substance]
        else:
            name = rd.substance_names[substance]
        ax.set_title('['+name+'] vs. t and x')

    return ax


def plot_bin_k_factors(rd, ax=None, indices=None):
    """
    Convenience function to inspect bin_k_factor in of
    ReactionDiffusion instance

    Parameters
    ----------
    rd: ReactionDiffusion
    ax: Axes instance or dict
        if ax is a dict it is used as **kwargs passed to
        matplotlib.pyplot.axes (default: None)
    indices: sequence of integers
        what factor sequences to plot
    """
    if not hasattr(ax, 'plot'):
        if ax is None:
            ax = plt.axes()
        else:
            ax = plt.axes(**ax)
    indices = indices or range(len(rd.bin_k_factor_span))
    factors = np.array(rd.bin_k_factor)
    x_edges = np.repeat(rd.x, 2)[1:]
    for i in indices:
        y_edges = np.pad(np.repeat(factors[:, i], 2), (0, 1), 'constant')
        ax.plot(x_edges, y_edges, label=i)
