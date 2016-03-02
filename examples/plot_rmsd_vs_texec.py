try:
    import cPickle as pickle
except ImportError:
    import pickle

import os
import struct

import numpy as np

from bokeh.io import curdoc
from bokeh.plotting import Figure, gridplot
from bokeh.models import ColumnDataSource, HoverTool

from rmsd_vs_texec import varied

results = pickle.load(open(
    os.path.join(os.path.dirname(__file__),
                 'analytic_diffusion_results.pkl'), 'rb'))

xkey, ykey = 'texec', 'rmsd_over_atol'
cols = [xkey, ykey, 'color'] + list(varied.keys())


hover = HoverTool(tooltips=[(k, '@'+k) for k in varied])

# Calculate colors
base_colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
               (1, 1, 0), (1, 0, 1), (0, 1, 1)]


def _color(params):
    hue_val, tone_val = params
    hue_vals, tone_vals = varied.values()
    tone = 1/4. + 3*tone_val/float(tone_vals[-1])/4.
    hue = base_colors[hue_vals.index(hue_val)]
    color = tuple(np.asarray(np.round(255*np.array(hue)*tone), dtype=int))
    return '#' + struct.pack('BBB', *color).encode('hex')

for k, v in results.items():
    results[k]['color'] = _color(k)


vals = results.values()
source = ColumnDataSource(data={k: [] for k in cols})
for k in cols:
    source.data[k] = [vals[idx][k] for idx in range(len(vals))]

top = Figure(plot_height=600, plot_width=800, title="%s vs. %s" % (xkey, ykey),
             x_axis_type="linear", y_axis_type="log",
             tools=[hover, 'pan', 'reset', 'box_zoom'])
top.xaxis.axis_label = xkey
top.yaxis.axis_label = ykey
cr = top.circle(x=xkey, y=ykey, source=source, size=10, color="color",
                line_color=None)


bottom = Figure(plot_height=600, plot_width=800, tools='pan,box_zoom')
bottom.line('tout', 'rmsd', source=source)


p = gridplot([[top], [bottom]])
curdoc().add_root(p)
