try:
    import cPickle as pickle
except ImportError:
    import pickle

import gzip
import os
import struct

import numpy as np

from bokeh.io import curdoc
from bokeh.plotting import Figure, gridplot
from bokeh.models import ColumnDataSource, HoverTool

from rmsd_vs_texec import varied


def main():
    results = pickle.load(gzip.open(
        os.path.join(os.path.dirname(__file__),
                     'analytic_diffusion_results.pkl'), 'rb'))

    xkey, ykey = 'texec', 'rmsd_over_atol'
    include_keys = list(varied.keys()) + [
        'nfev', 'njev',  'nprec_setup', 'nprec_solve', 'njacvec_dot',
        'nprec_solve_ilu', 'nprec_solve_lu',
        "n_steps",
        "n_rhs_evals",
        "n_lin_solv_setups",
        "n_err_test_fails",
        "n_nonlin_solv_iters",
        "n_nonlin_solv_conv_fails",
        "krylov_n_lin_iters",
        "krylov_n_prec_evals",
        "krylov_n_prec_solves",
        "krylov_n_conv_fails",
        "krylov_n_jac_times_evals",
        "krylov_n_iter_rhs"
    ]
    cols = [
        xkey, ykey, 'color',  # 'mrk'
    ] + include_keys

    hover = HoverTool(tooltips=[(k, '@'+k) for k in include_keys])

    # Calculate colors
    base_colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
                   (1, 1, 0), (1, 0, 1), (0, 1, 1)]

    def _color(params):
        hue_val, tone_val = params[:2]
        hue_vals, tone_vals = varied.values()[:2]
        tone = 1/4. + 3*tone_val/float(tone_vals[-1])/4.
        hue = base_colors[hue_vals.index(hue_val)]
        color = tuple(np.asarray(np.round(255*np.array(hue)*tone), dtype=int))
        return '#' + struct.pack('BBB', *color).encode('hex')

    for k, v in results.items():
        results[k]['color'] = _color(k)
        # results[k]['mrk'] = ('circle', 'square')[
        #     varied.values()[-1].index(k[-1])]

    vals = results.values()
    source = ColumnDataSource(data={k: [] for k in cols})
    for k in cols:
        source.data[k] = [vals[idx][k] for idx in range(len(vals))]

    top = Figure(plot_height=600, plot_width=800,
                 title="%s vs. %s" % (xkey, ykey),
                 x_axis_type="linear", y_axis_type="log",
                 tools=[hover, 'pan', 'reset', 'box_zoom'])
    top.xaxis.axis_label = xkey
    top.yaxis.axis_label = ykey
    top.scatter(
        x=xkey, y=ykey, source=source, size=10, color="color", line_color=None,
        # marker="mrk"
    )

    bottom = Figure(plot_height=600, plot_width=800, tools='pan,box_zoom')
    bottom.line('tout', 'rmsd', source=source)

    return gridplot([[top], [bottom]])


if __name__ == '__main__' or __name__.startswith('bk_'):
    curdoc().add_root(main())
