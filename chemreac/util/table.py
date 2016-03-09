# -*- coding: utf-8 -*-
"""
chemreac.util.table
-------------------

This module contains functions to generate tables.
"""
from __future__ import (absolute_import, division, print_function)

from chempy.units import get_derived_unit, to_unitless
from chempy.util.table import render_tex_to_pdf, tex_templates


def radyields2pdf_table(rd, output_dir=None, save=True, unit_registry=None,
                        siunitx=False, fmtstr='{0:.3f}', **kwargs):
    """ Generate a table with radiolytic yields

    Calls chempy.util.table.render_tex_to_pdf

    Parameters
    ----------
    rd: ReactionDiffusion
    output_dir: str
    save: bool
    unit_registry: dict
    siunitx: bool
    fmtstr: str
    \*\*kwargs:
        extends the table template dictionary
    """
    line_term = r' \\'
    col_delim = ' & '
    header = (col_delim.join(rd.substance_latex_names or rd.substance_names) +
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
    _pkgs = (['siunitx'] if siunitx else []) + [
        'booktabs', 'lscape', 'amsmath', 'hyperref']
    contents = tex_templates['document']['default'] % {
        'usepkg': '\n'.join([r'\usepackage{%s}' % pkg for pkg in _pkgs]),
        'begins': '\n'.join([r'\begin{%s}' % env for env in _envs]),
        'ends': '\n'.join([r'\end{%s}' % env for env in _envs[::-1]]),
        'table': table
    }
    return render_tex_to_pdf(contents, 'gvalues.tex', 'gvalues.pdf',
                             output_dir, save)
