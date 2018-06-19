import numpy as np
from .plotting import DEFAULT


def animate_C_vs_x(integr, ax, *, xunit=None, yunit=None, auto_ylim=False, title_fmt='t = {}',
                   ls=None, c=None, substances=None, yout=None):
    rd = integr.rd
    if yout is None:
        yout = integr.unitless_as('Cout', yunit)
    x_edges = np.repeat(integr.unitless_as('x', xunit), 2)[1:-1]
    y_edges = np.repeat(yout, 2, axis=1)

    if auto_ylim:
        ax.set_ylim([np.min(y_edges), np.max(y_edges)*1.05])

    ls = ls or DEFAULT['ls']
    c = c or DEFAULT['c']
    lines = []
    if rd.substance_latex_names in (None, [None]*rd.n):
        labels = rd.substance_names
    else:
        labels = ['$\\mathrm{'+n+'}$' for n in rd.substance_latex_names]

    for k in substances or rd.substance_names:
        si = rd.substance_names.index(k)
        lines.extend(ax.plot(x_edges, y_edges[0, :, si],
                             label=labels[si], ls=ls[si % len(ls)], c=c[si % len(c)]))
    tout_wu = integr.with_units('tout')
    ttl = ax.set_title(title_fmt.format(tout_wu[0]))

    def update(ti):
        ttl.set_text(title_fmt.format(tout_wu[ti]))
        for li, line in enumerate(lines):
            line.set_ydata(y_edges[ti, :, li])

    return update
