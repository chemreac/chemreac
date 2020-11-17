# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
import matplotlib.axes

from chemreac import ReactionDiffusion
from chemreac.integrate import run
from chemreac.util.plotting import (
    coloured_spy, plot_jacobian, plot_per_reaction_contribution,
    plot_C_vs_t_in_bin, plot_C_vs_x, plot_C_vs_t_and_x, plot_fields,
    plot_solver_linear_error, plot_solver_linear_excess_error,
)


def _get_decay_rd(N):
    return ReactionDiffusion(2, [[0]], [[1]], [3.14], N, D=[0.0, 0.0])


def _get_decay_Cref(N, y0, tout):
    y0 = np.asarray(y0)
    tout = tout.reshape((tout.size, 1))
    factor = np.exp(-(tout-tout[0])*3.14)
    Cref = np.hstack([
        np.hstack((y0[i*2]*factor, y0[i*2]*(1 - factor) + y0[i*2+1]))
        for i in range(N)])
    return Cref.reshape((tout.size, N, 2))


def test_get_decay_Cref():
    N = 3
    rd = _get_decay_rd(N)
    tout = np.linspace(0, 3.0, 7)
    y0 = [3.0, 1.0]*N
    integr = run(rd, y0, tout)
    Cref = _get_decay_Cref(N, y0, tout)
    assert np.allclose(Cref, integr.yout)


def test_plot_jacobian():
    # A -> B
    rd = _get_decay_rd(1)
    y0 = [3.0, 1.0]
    tout = np.linspace(0, 3.0, 7)
    integr = run(rd, y0, tout)
    axes = plot_jacobian(rd, integr.tout, integr.yout, [0, 1])
    for ax in axes:
        assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_per_reaction_contribution():
    rd = _get_decay_rd(1)
    tout = np.linspace(0, 3.0, 7)
    y0 = [3.0, 1.0]
    integr = run(rd, y0, tout)
    axes = plot_per_reaction_contribution(integr, [0, 1])
    for ax in axes:
        assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_C_vs_t_in_bin():
    N = 3
    rd = _get_decay_rd(N)
    tout = np.linspace(0, 3.0, 7)
    y0 = [3.0, 1.0]*N
    integr = run(rd, y0, tout)
    ax = plot_C_vs_t_in_bin(rd, tout, integr.yout, 0)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_C_vs_x():
    N = 3
    rd = _get_decay_rd(N)
    tout = np.linspace(0, 3.0, 7)
    y0 = [3.0, 1.0]*N
    integr = run(rd, y0, tout)
    ax = plot_C_vs_x(rd, tout, integr.yout, [0, 1], 6)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_C_vs_t_and_x():
    import mpl_toolkits
    N = 3
    rd = _get_decay_rd(N)
    tout = np.linspace(0, 3.0, 7)
    y0 = [3.0, 1.0]*N
    integr = run(rd, y0, tout)
    ax = plot_C_vs_t_and_x(rd, tout, integr.yout, 0)
    assert isinstance(ax, mpl_toolkits.mplot3d.Axes3D)


def test_plot_fields():
    # modulation in x means x_center
    # A -> B # mod1 (x**2)
    # C -> D # mod1 (x**2)
    # E -> F # mod2 (sqrt(x))
    # G -> H # no modulation
    k = np.array([3.0, 7.0, 13.0, 22.0])
    N = 5
    n = 8
    D = np.zeros(n)
    x = np.linspace(3, 7, N+1)
    xc = x[:-1] + np.diff(x)/2
    fields = [
        [xc[i]*xc[i] for i in range(N)],
        [xc[i]*xc[i] for i in range(N)],
        [xc[i]**0.5 for i in range(N)]
    ]
    g_values = [
        [-k[0], k[0], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -k[1], k[1], 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -k[2], k[2], 0.0, 0.0],
    ]
    rd = ReactionDiffusion(
        n,
        [[6]],
        [[7]],
        k=k[3:], N=N, D=D, x=x, fields=fields,
        g_values=g_values
    )
    ax = plot_fields(rd, indices=[0, 2])
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_solver_linear_error():
    N = 3
    rd = _get_decay_rd(N)
    tout = np.linspace(0, 3.0, 7)
    y0 = [3.0, 1.0]*N
    integr = run(rd, y0, tout)
    Cref = _get_decay_Cref(N, y0, tout)
    ax = plot_solver_linear_error(integr, Cref)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_plot_solver_linear_excess_error():
    N = 3
    rd = _get_decay_rd(N)
    tout = np.linspace(0, 3.0, 7)
    y0 = [3.0, 1.0]*N
    integr = run(rd, y0, tout)
    Cref = _get_decay_Cref(N, y0, tout)
    ax = plot_solver_linear_excess_error(integr, Cref)
    assert isinstance(ax, matplotlib.axes.Axes)


def test_coloured_spy():
    from matplotlib.axes import Axes
    A = np.arange(9).reshape((3, 3))
    for log in (False, True, -5):
        ax_im, ax_cb = coloured_spy(A, log=log)
        assert isinstance(ax_im, Axes)
        assert isinstance(ax_cb, Axes)
