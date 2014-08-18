# -*- coding: utf-8 -*-

"""
Convenince functions for representing reactin systems in tables.
"""

import subprocess
import tempfile

def rsys2tablines(rsys, substances, rref0=1, coldelim=' & ',
                  tex=True, rxnarrow=r'$\rightarrow$', ref_fmt='{}'):
    param_fmt = '{0:.3g}'  # Could be taken from instance of Reaction (subclass)
    _get_name = lambda sn: getattr(substances[sn], 'tex_name' if tex else 'name')
    lines = []
    for ri, rxn in enumerate(rsys.rxns):
        lines.append(coldelim.join([
            str(rref0+ri),
            ' + '.join([('' if num == 1 else str(num)) + _get_name(sn) for
                        sn, num in rxn.reactants.items()]),
            rxnarrow,
            ' + '.join([('' if num == 1 else str(num)) + _get_name(sn) for
                       sn, num in rxn.products.items()]),
            param_fmt.format(rxn.k),
            ref_fmt.format(rxn.ref)
        ]))
    return lines


def rsys2table(rsys, substances, table_template=None, table_template_dict=None, **kwargs):
    """
    Renders user provided table_template with table_template_dict which also has 'body'
    entry generated from `rsys2tablines`.

    Defaults is LaTeX table requiring booktabs package to be used
    (add \usepackage{booktabs} to preamble).
    """
    line_term = r' \\'
    if table_template is None:
        table_template = r"""
\begin{%(table_env)s}
\centering
\label{tab:%(label)s}
\caption[%(short_cap)s]{%(long_cap)s}
\begin{tabular}{%(alignment)s}
\toprule
%(header)s
\midrule
%(body)s
\bottomrule
\end{tabular}
\end{%(table_env)s}"""

    defaults = {
        'table_env': 'table',  # e.g. longtable
        'alignment': 'llllll',  # e.g. llllSl for siunitx
        'header': kwargs.get('coldelim', ' & ').join([
            'Id.', 'Reactants', '','Products', 'Rate constant', 'Ref'
        ]) + line_term,
        'short_cap': rsys.name,
        'long_cap': rsys.name,
        'label': (rsys.name or 'None').lower()
    }

    if table_template_dict is None:
        table_template_dict = defaults
    else:
        for k, v in defaults:
            if not k in table_template_dict:
                table_template_dict[k] = v

    if 'body' in table_template_dict:
        raise KeyError("There is already a 'body' key in table_template_dict")
    table_template_dict['body'] = (line_term + '\n').join(rsys2tablines(
        rsys, substances, **kwargs)) + line_term
    return table_template % table_template_dict


def rsys2pdf_table(rsys, substances, output_dir=None, tex_template=None,
             tex_template_dict=None, **kwargs):
    """
    kwargs passed on to `rsys2table`
    """
    if tex_template is None:
        tex_template = r"""
\documentclass{article}
\pagestyle{empty}
\usepackage{booktabs}
\usepackage{lscape}
\begin{document}
%(begins)s
%(table)s
%(ends)s
\end{document}
"""

    envs = ['landscape', 'tiny']
    defaults = {
        'begins': '\n'.join([r'\begin{%s}' % env for env in envs]),
        'ends': '\n'.join([r'\end{%s}' % env for env in envs[::-1]])
    }

    if tex_template_dict is None:
        tex_template_dict = defaults
    else:
        for k, v in defaults:
            if not k in tex_template_dict:
                tex_template_dict[k] = v

    if 'table' in tex_template_dict:
        raise KeyError("There is already a 'tex' key in tex_template_dict")
    tex_template_dict['table'] = rsys2table(rsys, substances, **kwargs)

    contents = tex_template % tex_template_dict

    with tempfile.NamedTemporaryFile('wt', suffix='.tex', dir=output_dir) as tmpfh:
        cmds = ['pdflatex', tmpfh.name]
        tmpfh.write(contents)
        tmpfh.flush()
        p = subprocess.Popen(cmds + [tmpfh.name], cwd=output_dir)
        retcode = p.wait()
        if retcode:
            raise RuntimeError(' '.join(cmds) + "\n returned with exit status {}".format(
                retcode))
