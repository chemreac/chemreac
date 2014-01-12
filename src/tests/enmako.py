#!/usr/bin/env python

from mako.template import Template
from mako.exceptions import text_error_template
import argh

def render_mako_template_to(template_path, outpath, subsd):
    template_str = open(template_path, 'rt').read()
    with open(outpath, 'wt') as ofh:
        try:
            rendered = Template(template_str, input_encoding='utf-8',
                                output_encoding='utf-8').render(**subsd)
        except:
            print(text_error_template().render())
            raise
        ofh.write(rendered)

def enmako(template_path, gen_subsd_eval=None, shell_cmd_subs=None, json_subs=None, pickle_subs=None, outpath=None):
    """
    User provides a template file path, e.g. `index.html.mako`
    and either a python expression (gen_subsd_eval) which evaluates to a dict e.g. '{"title": "Welcome"}'
    or a shell command (shell_cmd_subs) which outputs lines of form: title=Welcome
    If outpath is not given, template_path will be stripped from trailing `.mako` (required in that case)

    Note: json does not support integer keys in dicts
    """
    if outpath == None:
        assert template_path.endswith('.mako')
        outpath = template_path[:-5]
    assert sum([0 if x==None else 1 for x in [gen_subsd_eval, shell_cmd_subs, json_subs, pickle_subs]]) == 1
    if gen_subsd_eval:
        subsd = eval(gen_subsd_eval)
    if json_subs:
        import json
        subsd = json.load(open(json_subs, 'rt'))
    if pickle_subs:
        try:
            import cPickle as pickle
        except ImportError:
            import pickle
        subsd = pickle.load(open(pickle_subs, 'rt'))
    if shell_cmd_subs:
        import subprocess
        outp = subprocess.check_output(shell_cmd_subs.split())
        subsd = dict([x.split('=') for x in outp.split('\n')[:-1]])
    render_mako_template_to(template_path, outpath, subsd)

if __name__ == '__main__':
    argh.dispatch_command(enmako)
