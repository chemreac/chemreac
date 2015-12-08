# -*- coding: utf-8 -*-

"""
chemreac.util.table
-------------------

Convenience functions for presenting reaction systems in tables.

"""

import os
import shutil
import subprocess
import tempfile

from chemreac.units import to_unitless, get_derived_unit

tex_templates = {
    'document': {
        'default': r"""
\documentclass{article}
\pagestyle{empty}
%(usepkg)s
\begin{document}
%(begins)s
%(table)s
%(ends)s
\end{document}
"""
    },
    'table': {
        'default': r"""
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
\end{%(table_env)s}""",
        'longtable': r"""
\begin{%(table_env)s}{%(alignment)s}
\caption[%(short_cap)s]{%(long_cap)s
\label{tab:%(label)s}}\\
\toprule
%(header)s
\midrule
%(body)s
\bottomrule
\end{%(table_env)s}"""
    }
}


def render_tex_to_pdf(contents, texfname, pdffname, output_dir, save):
    created_tempdir = False
    try:
        if output_dir is None:
            output_dir = tempfile.mkdtemp()
            created_tempdir = True
        texpath = os.path.join(output_dir, texfname)
        pdfpath = os.path.join(output_dir, pdffname)
        cmds = ['pdflatex', '-halt-on-error', '-interaction',
                'batchmode', texfname]
        with open(texpath, 'wt') as ofh:
            ofh.write(contents)
            ofh.flush()
        with open(pdfpath + '.out', 'wb') as logfile:
            p = subprocess.Popen(cmds, cwd=output_dir,
                                 stdout=logfile, stderr=logfile)
            retcode = p.wait()
            p = subprocess.Popen(cmds, cwd=output_dir,
                                 stdout=logfile, stderr=logfile)
            retcode += p.wait()
        if retcode:
            fmtstr = "{}\n returned with exit status {}"
            raise RuntimeError(fmtstr.format(' '.join(cmds), retcode))
        else:
            return pdfpath
    finally:
        if save is True or save == 'True':
            pass
        else:
            if save is False or save == 'False':
                if created_tempdir:
                    shutil.rmtree(output_dir)
            else:
                # interpret path to copy pdf to.
                shutil.copy(pdfpath, save)


def rsys2tablines(rsys, substances, rref0=1, coldelim=' & ',
                  tex=True, rxnarrow=r'$\rightarrow$', ref_fmt='{}',
                  unit_registry=None, unit_fmt='{}'):
    """
    Generates a table representation of a ReactionSystem.

    Parameters
    ----------
    rsys: ReactionSystem
    substances: sequence of Substance instances
    rref0: integer
        default start of index counter (default: 1)
    coldelim: string
        column delimiter (default: ' & ')
    tex: bool
        use latex formated output (default: True)
    rxnarrow: string
        default: '\$\\rightarrow\$'
    ref_fmt: string or callable
        format string of ``ref`` attribute of reactions
    unit_registry: unit registry
        optional (default: None)
    """

    param_fmt = '{0:.3g}'  # Could be taken from Reaction instance

    def _get_name(sn):
        return getattr(substances[sn], 'tex_name' if tex else 'name')
    lines = []
    for ri, rxn in enumerate(rsys.rxns):
        if unit_registry is not None:
            kunit = (get_derived_unit(unit_registry,
                                      'concentration')**(1-rxn.order) /
                     get_derived_unit(unit_registry, 'time'))
            k = to_unitless(rxn.k, kunit)
        else:
            kunit = 1
            k = rxn.k
        lines.append(coldelim.join([
            str(rref0+ri),
            ' + '.join([('' if num == 1 else str(num)) + _get_name(sn) for
                        sn, num in rxn.reactants.items() if num > 0]),
            rxnarrow,
            ' + '.join([('' if num == 1 else str(num)) + _get_name(sn) for
                        sn, num in rxn.products.items() if num > 0]),
            param_fmt.format(k),
            unit_fmt.format(kunit),
            ref_fmt(rxn.ref) if callable(ref_fmt) else ref_fmt.format(rxn.ref)
        ]))
    return lines


def rsys2table(rsys, substances, table_template=None,
               table_template_dict=None, **kwargs):
    r"""
    Renders user provided table_template with table_template_dict which
    also has 'body' entry generated from `rsys2tablines`.

    Defaults is LaTeX table requiring booktabs package to be used
    (add \usepackage{booktabs} to preamble).

    Parameters
    ==========
    rsys: ReactionSystem
    substances: sequence of Substance instances
    table_template: string
    table_tempalte_dict: dict used to render table_template (excl. "body")
    longtable: bool
        use longtable in defaults. (default: False)
    **kwargs:
        passed onto rsys2tablines
    """
    siunitx = kwargs.pop('siunitx', False)
    line_term = r' \\'
    defaults = {
        'table_env': 'longtable' if kwargs.pop(
            'longtable', False) else 'table',
        'alignment': 'llllSll' if siunitx else 'lllllll',
        'header': kwargs.get('coldelim', ' & ').join([
            'Id.', 'Reactants', '', 'Products', '{Rate constant}',
            'Unit', 'Ref'
        ]) + line_term,
        'short_cap': rsys.name,
        'long_cap': rsys.name,
        'label': (rsys.name or 'None').lower()
    }

    if table_template_dict is None:
        table_template_dict = defaults
    else:
        for k, v in defaults:
            if k not in table_template_dict:
                table_template_dict[k] = v

    if 'body' in table_template_dict:
        raise KeyError("There is already a 'body' key in table_template_dict")
    table_template_dict['body'] = (line_term + '\n').join(rsys2tablines(
        rsys, substances, **kwargs)) + line_term

    if table_template is None:
        if table_template_dict['table_env'] == 'longtable':
            table_template = tex_templates['table']['longtable']
        else:
            table_template = tex_templates['table']['default']

    return table_template % table_template_dict


def rsys2pdf_table(rsys, substances, output_dir=None, doc_template=None,
                   doc_template_dict=None, save=True, landscape=False,
                   **kwargs):
    """
    Convenience function to render a ReactionSystem as
    e.g. a pdf using e.g. pdflatex.

    Parameters
    ==========
    rsys: ReactionSystem
    substances: sequence of Substance instances
    output_dir: path to output directory
        (default: system's temporary folder)
    doc_template: string
        LaTeX boiler plate temlpate including preamble,
        document environment etc.
    doc_template_dict: dict (string -> string)
        dict used to render temlpate (excl. 'table')
    longtable: bool
        use longtable in defaults. (default: False)
    **kwargs:
        passed on to `rsys2table`
    """
    if doc_template is None:
        doc_template = tex_templates['document']['default']
    _pkgs = ['booktabs'] + (['lscape'] if landscape else [])
    if kwargs.get('longtable', False):
        _pkgs += ['longtable']
    if kwargs.get('siunitx', False):
        _pkgs += ['siunitx']
    _envs = ['tiny'] + (['landscape'] if landscape else [])
    defaults = {
        'usepkg': '\n'.join([r'\usepackage{%s}' % pkg for pkg in _pkgs]),
        'begins': '\n'.join([r'\begin{%s}' % env for env in _envs]),
        'ends': '\n'.join([r'\end{%s}' % env for env in _envs[::-1]])
    }

    if doc_template_dict is None:
        doc_template_dict = defaults
    else:
        for k, v in defaults:
            if k not in doc_template_dict:
                doc_template_dict[k] = v

    if 'table' in doc_template_dict:
        raise KeyError("There is already a 'table' key in doc_template_dict")
    doc_template_dict['table'] = rsys2table(rsys, substances, **kwargs)

    contents = doc_template % doc_template_dict

    texfname = (rsys.name or 'output') + '.tex'
    pdffname = (rsys.name or 'output') + '.pdf'
    return render_tex_to_pdf(contents, texfname, pdffname, output_dir, save)


def radyields2pdf_table(rd, output_dir=None, save=True, unit_registry=None,
                        siunitx=False, fmtstr='{0:.3f}', **kwargs):
    line_term = r' \\'
    col_delim = ' & '
    header = (col_delim.join(rd.substance_tex_names or rd.substance_names) +
              line_term)
    lines = []
    for cur_gs in rd.g_values:
        if unit_registry is not None:
            gunit = get_derived_unit(unit_registry, 'radiolytic_yield')
            cur_gs = to_unitless(cur_gs, gunit)
        lines.append(col_delim.join(map(
            lambda v: fmtstr.format(v), cur_gs)) + line_term)
    table_template_dict = {
        'table_env': 'table',
        'alignment': ('@{}S' if siunitx else '@{}l')*rd.n,
        'header': header,
        'short_cap': 'G-values',
        'long_cap': 'G-values',
        'label': 'none',
        'body': '\n'.join(lines)
    }
    table_template_dict.update(kwargs)
    table = tex_templates['table']['default'] % table_template_dict

    _envs = ['landscape', 'tiny']
    _pkgs = (['siunitx'] if siunitx else []) + ['booktabs', 'lscape']
    contents = tex_templates['document']['default'] % {
        'usepkg': '\n'.join([r'\usepackage{%s}' % pkg for pkg in _pkgs]),
        'begins': '\n'.join([r'\begin{%s}' % env for env in _envs]),
        'ends': '\n'.join([r'\end{%s}' % env for env in _envs[::-1]]),
        'table': table
    }
    return render_tex_to_pdf(contents, 'gvalues.tex', 'gvalues.pdf',
                             output_dir, save)
