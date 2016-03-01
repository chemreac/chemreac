try:
    import cPickle as pickle
except ImportError:
    import pickle

from bokeh.io import curdoc
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HoverTool

from rmsd_vs_texec import varied

results = pickle.load(open('analytic_diffusion_results.pkl', 'rt'))

xkey, ykey = 'texec', 'rmsd_over_atol'
cols = [xkey, ykey] + list(varied.keys())

source = ColumnDataSource(data={k: [] for k in cols})

hover = HoverTool(tooltips=[(k, '@'+k) for k in varied])

vals = results.values()
for k in cols:
    source.data[k] = [vals[idx][k] for idx in range(len(vals))]

p = Figure(plot_height=600, plot_width=800, title="", toolbar_location=None, tools=[hover, 'pan'])
p.circle(x=xkey, y=ykey, source=source, size=7, line_color=None) # color="color",

#right = Figure(x='t', y='y')

curdoc().add_root(p)
