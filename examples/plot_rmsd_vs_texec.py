from __future__ import (absolute_import, division, print_function)

try:
    import cPickle as pickle
except ImportError:
    import pickle

import gzip
import os
import struct

import numpy as np

xkey, ykey = 'times_cpu', 'rmsd_over_atol'


def read():
    basename = os.path.splitext(os.path.basename(__file__))[0]
    assert basename.startswith('plot_')
    basename = basename[len('plot_'):]
    source = os.path.join(os.path.dirname(__file__), basename + '.pkl')
    if not os.path.exists(source):
        raise IOError("%s does not exist. Run rmsd_vs_texec.py first" % source)
    results = pickle.load(gzip.open(source, 'rb'))
    varied_keys = results.pop('varied_keys')
    varied_vals = results.pop('varied_values')
    if len(varied_keys) != 3 or len(varied_vals) != 3:
        raise ValueError("Script assumes 3 parameters (hue, tone, marker)")

    base_colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
                   (1, 1, 0), (1, 0, 1), (0, 1, 1)]

    def _color(params):
        hue_val, tone_val = params[:2]
        hue_vals, tone_vals = varied_vals[:2]
        tone = 1/4. + 3*tone_val/float(tone_vals[-1])/4.
        hue = base_colors[hue_vals.index(hue_val)]
        color = tuple(np.asarray(np.round(255*np.array(hue)*tone), dtype=int))
        return '#' + struct.pack('BBB', *color).encode('hex')

    for k, v in results.items():
        results[k]['color'] = _color(k)

    return results, varied_keys, varied_vals


def get_bokeh_fig():
    from bokeh.plotting import Figure  # , gridplot
    from bokeh.models import ColumnDataSource, HoverTool
    results, varied_keys, varied_vals = read()
    include_keys = varied_keys + [
        'nfev', 'njev',  'nprec_setup', 'nprec_solve', 'njacvec_dot',
        'nprec_solve_ilu', 'nprec_solve_lu', "n_steps", "n_rhs_evals",
        "n_lin_solv_setups", "n_err_test_fails", "n_nonlin_solv_iters",
        "n_nonlin_solv_conv_fails", "krylov_n_lin_iters",
        "krylov_n_prec_evals", "krylov_n_prec_solves", "krylov_n_conv_fails",
        "krylov_n_jac_times_evals", "krylov_n_iter_rhs"
    ]
    cols = [xkey, ykey, 'color'] + include_keys
    sources = {}
    varied3 = varied_vals[2]
    keys = list(results.keys())
    vals = list(results.values())
    for val in varied3:
        sources[val] = ColumnDataSource(data={k: [] for k in cols})
        for k in cols:
            sources[val].data[k] = [vals[idx].get(k, None) for idx in
                                    range(len(vals)) if keys[idx][2] == val]
    hover = HoverTool(tooltips=[(k, '@'+k) for k in include_keys])
    top = Figure(
        plot_height=600, plot_width=800, title="%s vs. %s" % (ykey, xkey),
        x_axis_type="linear", y_axis_type="log", tools=[
            hover, 'pan', 'reset', 'box_zoom', 'wheel_zoom', 'save'])
    top.xaxis.axis_label = xkey
    top.yaxis.axis_label = ykey
    for source, marker in zip(sources.values(), ['circle', 'diamond']):
        top.scatter(x=xkey, y=ykey, source=source, size=9, color="color",
                    line_color=None, marker=marker)
    return top


def plot_with_matplotlib(savefig='none', dpi=300, errorbar=False, linx=False,
                         liny=False):
    import matplotlib.pyplot as plt
    results, varied_keys, varied_vals = read()
    plt.rc('font', family='serif')

    def label(data):
        if data[varied_keys[1]] == varied_vals[1][-1]:
            meth = data['method']
            meth = meth.replace('bdf', r'$\mathrm{BDF}$')
            meth = meth.replace('adams', r'$\mathrm{Adams}$')
            return '$%d$, %s' % (data['nstencil'], meth)

    for params, data in results.items():
        mrkr = 'od'[varied_vals[2].index(params[2])]
        if errorbar:
            cb, kwargs = plt.errorbar, dict(xerr=2*data['d'+xkey])
        else:
            cb, kwargs = plt.plot, {}
        cb(data[xkey], data[ykey], marker=mrkr, color=data['color'],
           label=label(data), ls='None', **kwargs)
    ax = plt.gca()
    if not linx:
        ax.set_xscale('log')
    if not liny:
        ax.set_yscale('log')
    plt.xlabel('$t_{exec}$')
    plt.ylabel(r'$\mathrm{RMSD}\ /\ \mathrm{atol}$')

    handles, labels = ax.get_legend_handles_labels()

    # reverse the order
    ax.legend(handles[::-1], labels[::-1])

    # or sort them by labels
    import operator
    hl = sorted(zip(handles, labels),
                key=operator.itemgetter(1))
    handles2, labels2 = zip(*hl)

    ax.legend(handles2, labels2, numpoints=1, prop={'size': 12})
    if savefig.lower() == 'none':
        plt.show()
    else:
        plt.savefig(savefig, dpi=dpi)

if __name__.startswith('bk_'):
    # e.g:
    #
    #  $ bokeh html plot_rmsd_vs_texec.py
    #  $ bokeh serve plot_rmsd_vs_texec.py
    from bokeh.io import curdoc
    curdoc().add_root(get_bokeh_fig())

if __name__ == '__main__':
    # e.g:
    #
    #  $ python plot_rmsd_vs_texec.py
    #  $ python plot_rmsd_vs_texec.py --savefig figure.png --dpi 100
    import argh
    argh.dispatch_command(plot_with_matplotlib)
