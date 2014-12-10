# -*- coding: utf-8 -*-

"""
plotting
--------

convenience functions to create matplotlib plots
of results.
"""

import numpy as np
from chemreac.util.analysis import solver_linear_error_from_integration
from chemreac.util.banded import get_jac_row_from_banded
from chemreac.util.pyutil import dict_with_defaults
import matplotlib.pyplot as plt


def _init_axes(ax=None):
    if ax is None:
        ax = plt.axes()
    else:
        if not hasattr(ax, 'plot'):
            ax = plt.axes(**ax)
    return ax


def save_and_or_show_plot(show=None, savefig='None'):
    """
    Convenience method for either showing or saving current matplotlib
    figure.

    Parameters
    ----------
    show: bool or None
        Show plot, when None only show when savefig is not used
        default: None
    savefig: string
        path to output file of figure. If extension is html, mpld3
        will be used to generate a d3 backed html output.

    """
    if savefig is not None and savefig != 'None':
        if savefig.endswith('.html'):
            # Export using mpld3
            import mpld3
            open(savefig, 'wt').write(mpld3.fig_to_html(plt.gcf()))
        else:
            plt.savefig(savefig)

        if show is None:
            show = False
    else:
        if show is None:
            show = True

    if show:
        plt.show()


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
    legend_kwargs = dict_with_defaults(legend_kwargs,
                                       dict(loc=None, prop={'size': 11}))
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
    print_names = rd.substance_tex_names or rd.substance_names
    axes = _plot_analysis(_get_jac_row_over_t, print_names, rd, tout, yout,
                          indices, titles=[
                              print_names[i] for i in indices], **kwargs)
    for ax in axes:
        ax.set_ylabel(
            "$\\frac{\\partial r_{tot}}{\partial C_i}~/~s^{-1}$")
    return axes


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
    print_names = rd.substance_tex_names or rd.substance_names
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
        substance_idxs = range(rd.n)
    else:
        substance_idxs = [
            s if isinstance(s, int) else rd.substance_names.index(s) for
            s in substances]

    if labels is None:
        try:
            names = rd.substance_tex_names or rd.substance_names
        except AttributeError:
            names = list(map(str, substance_idxs))
        labels = [names[i] for i in substance_idxs]
    else:
        assert len(labels) == len(substance_idxs)
    return ax, substance_idxs, labels


def plot_C_vs_t_in_bin(
        rd, tout, yout, bi=0, ax=None, labels=None,
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
    legend_kwargs = dict_with_defaults(legend_kwargs,
                                       dict(loc='best', prop={'size': 11}))
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
    """
    Plots concentration as function of x for selected
    substances at time index 'ti'.

    Parameters
    ----------
    rd: ReactionDiffusion
    tout: 1D array of floats
    yout: output from solver
    substances: sequence of indies or names of substances
    ti: int
        time index
    ax: Axes instance
    labels: sequence of strings
    xscale: matplotlib scale choice (e.g. 'log', 'symlog')
    yscale: matplotlib scale choice (e.g. 'log', 'symlog')
    basetitle: string

    Returns
    =======
    Axes instance
    """
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
    return ax


def plot_C_vs_t_and_x(rd, tout, yout, substance, ax=None, log10=False,
                      **plot_kwargs):
    """
    Plots 3D surface of concentration as function of time and x for a
    selected substance.

    Parameters
    ----------
    rd: ReactionDiffusion
    tout: 1D array of floats
    yout: output from solver
    substance: int or string
        index or name of substance
    ax: Axes instance
    log10: bool
        Use log logarithmic (base 10) axis
    **plot_kwargs:
        passed onto plot_surface

    Returns
    =======
    Axes3D instance
    """
    # it would be nice to accpet kwargs
    #    xscale='log', yscale='log', zscale='log'
    # but it's currently not supported by matplotlib:
    # http://matplotlib.1069221.n5.nabble.com/
    #     plot-surface-fails-with-log-axes-td10206.html
    substance = (substance if isinstance(substance, int) else
                 rd.substance_names.index(substance))
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    ax = ax or plt.subplot(1, 1, 1, projection='3d')
    assert isinstance(ax, Axes3D)

    xty = [rd.xcenters, tout, yout]
    x_, t_, y_ = list(map(np.log10, xty)) if log10 else xty
    X, T = np.meshgrid(x_, t_)
    if 'cmap' not in plot_kwargs:
        plot_kwargs['cmap'] = cm.gist_earth
    ax.plot_surface(X, T, y_[:, :, substance],
                    **plot_kwargs)

    fmtstr = "$log_{{10}}$({})" if log10 else "{}"
    ax.set_xlabel(fmtstr.format('x / m'))
    ax.set_ylabel(fmtstr.format('time / s'))
    ax.set_zlabel(fmtstr.format('C / M'))
    if rd.substance_names:
        if rd.substance_tex_names:
            name = rd.substance_tex_names[substance]
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
        if ax is a dict it is used as keyword arguments passed to
        matplotlib.pyplot.axes (default: None)
    indices: sequence of integers
        what factor sequences to plot
    """
    ax = _init_axes(ax)
    indices = indices or range(len(rd.bin_k_factor_span))
    factors = np.array(rd.bin_k_factor)
    x_edges = np.repeat(rd.x, 2)[1:]
    for i in indices:
        y_edges = np.pad(np.repeat(factors[:, i], 2), (0, 1), 'constant')
        ax.plot(x_edges, y_edges, label=i)
    return ax


def plot_solver_linear_error(
        integration, Cref=0, ax=None, x=None, ti=slice(None), bi=0, si=0,
        plot_kwargs=None, fill_between_kwargs=None, scale_err=1.0, fill=True,
        **kwargs):
    """
    Parameters
    ----------
    integration: chemreac.integrate.Integration
        result from integration.
    Cref: array or float
        analytic solution to compare with
    ax: Axes instance or dict
        if ax is a dict it is used as key word arguments passed to
        matplotlib.pyplot.axes (default: None)
    x: array
        (optional) x-values, when None it is deduced to be
        either t or x (when ti or bi are slices repecitvely)
        (default: None)
    ti: slice
        time indices
    bi: slice
        bin indices
    si: integer
        specie index
    plot_kwargs: dict
        keyword arguments passed to matplotlib.pyplot.plot (default: None)
    fill_between_kwargs: dict
        keyword arguments passed to matplotlib.pyplot.fill_between
        (default: None)
    scale_err: float
        value with which errors are scaled. (default: 1.0)
    fill: bool
        whether or not to fill error span
    \*\*kwargs
        common keyword arguments of plot_kwargs and fill_between_kwargs,
        e.g. 'color', (default: None).

    See Also
    --------
    plot_solver_linear_excess_error
    """
    ax = _init_axes(ax)
    Cerr = integration.Cout - Cref
    if x is None:
        if isinstance(ti, slice):
            x = integration.tout[ti]
        elif isinstance(bi, slice):
            x = integration.rd.xcenters[bi]
        else:
            raise NotImplementedError("Failed to deduce x-axis.")

    plt.plot(x, scale_err*Cerr[ti, bi, si], **dict_with_defaults(
        plot_kwargs, kwargs))

    if fill:
        le_l, le_u = solver_linear_error_from_integration(
            integration, ti=ti, bi=bi, si=si)
        Cerr_u = le_u - Cref[ti, bi, si]
        Cerr_l = le_l - Cref[ti, bi, si]
        plt.fill_between(x, scale_err*Cerr_l,
                         scale_err*Cerr_u, **dict_with_defaults(
                             fill_between_kwargs, {'alpha': 0.2}, kwargs))
    return ax


def plot_solver_linear_excess_error(integration, Cref, ax=None, x=None,
                                    ti=slice(None), bi=0, si=0, **kwargs):
    """
    Plots the excess error commited by the intergrator, divided by the span
    of the tolerances (atol + rtol*|y_i|).

    Parameters
    ----------
    integration: chemreac.integrate.Integration
        result from integration.
    Cref: array or float
        analytic solution to compare with
    ax: Axes instance or dict
        if ax is a dict it is used as \*\*kwargs passed to
        matplotlib.pyplot.axes (default: None)
    x: array
        (optional) x-values, when None it is deduced to be
        either t or x (when ti or bi are slices repecitvely)
        (default: None)
    ti: slice
        time indices
    bi: slice
        bin indices
    si: integer
        specie index
    plot_kwargs: dict
        keyword arguments passed to matplotlib.pyplot.plot (default: None)
    fill_between_kwargs: dict
        keyword arguments passed to matplotlib.pyplot.fill_between
        (default: None)
    scale_err: float
        value with which errors are scaled. (default: 1.0)
    fill: bool
        whether or not to fill error span
    \*\*kwargs:
        common keyword arguments of plot_kwargs and fill_between_kwargs,
        e.g. 'color', (default: None).

    See Also
    --------
    plot_solver_linear_error
    """
    if x is None:
        if isinstance(ti, slice):
            x = integration.tout[ti]
        elif isinstance(bi, slice):
            x = integration.rd.xcenters[bi]
        else:
            raise NotImplementedError("Failed to deduce x-axis.")
    ax = _init_axes(ax)
    le_l, le_u = solver_linear_error_from_integration(integration, ti, bi, si)
    Eexcess_l = Cref[ti, bi, si] - le_l  # Excessive if negative
    Eexcess_u = Cref[ti, bi, si] - le_u  # Excessive if positive
    Eexcess_l[np.argwhere(Eexcess_l >= 0)] = 0
    Eexcess_u[np.argwhere(Eexcess_u <= 0)] = 0
    fused = np.concatenate((Eexcess_l[..., np.newaxis],
                            Eexcess_u[..., np.newaxis]), axis=-1)
    indices = np.argmax(abs(fused), axis=-1)
    Eexcess = fused[np.indices(indices.shape), indices][0, ...]
    le_span = le_u - le_l
    ax.plot(integration.tout, Eexcess/le_span, **kwargs)
    return ax
