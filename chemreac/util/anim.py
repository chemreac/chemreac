import numpy as np
from matplotlib import animation as manim

def animate_C_vs_x(integr, ax, xunit=None, yunit=None, auto_ylim=False, title_fmt='t = {}'):
    rd = integr.rd
    xscal = to_unitless(get_derived_unit(rd.unit_registry, 'length'), xunit) if xunit else 1
    yscal = to_unitless(get_derived_unit(rd.unit_registry, 'concentration'), xunit) if yunit else 1
    if auto_ylim:
        ylim = ax.get_ylim()
        ax.set_ylim([np.min(integr.Cout)*yscal, np.max(integr.Cout)*yscal*1.05])

    lines = []

    if rd.substance_latex_names in (None, [None]*rd.n):
        labels = rd.substance_names
    else:
        labels = ['$'+n+'$' for n in rd.substance_latex_names]

    x_edges = np.repeat(rd.x, 2)[1:-1]
    y_edges = np.repeat(integr.Cout, 2, axis=1)
    for si, k in enumerate(rd.substance_names):
        lines.extend(ax.plot(x_edges*xscal, y_edges[0, :, si]*yscal, label=labels[si]))
    tout_wu = integr.with_units('tout')
    ttl = ax.set_title(title_fmt.format(tout_wu[0]))
    def update(ti):
        ttl.set_text(title_fmt.format(tout_wu[ti]))
        for li, line in enumerate(lines):
            line.set_ydata(y_edges[ti, :, li]*yscal)
        #return lines

    return update
