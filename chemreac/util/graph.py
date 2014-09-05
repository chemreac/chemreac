# -*- coding: utf-8 -*-

import subprocess
import tempfile

"""
Convenince functions for representing reaction systems as graphs.
"""


def rsys2dot(rsys, substances, tex=False, rprefix='r', rref0=1,
             nodeparams='[label={} shape=diamond]'):
    """
    Returns list of lines of DOT (graph description language)
    formated graph.

    Parameters
    ==========
    rsys: ReactionSystem
    substances: sequence of Substance instances
    tex: bool (default False)
        If set True, output will be LaTeX formated
    (Substance need to have tex_name attribute set)
    rprefix: string
        Reaction enumeration prefix, default: r
    rref0: integer
        Reaction enumeration inital counter value, default: 1
    nodeparams: string
        DOT formated param list, default: [label={} shape=diamond]

    Returns
    =======
    list of lines of DOT representation of the graph representation.

    """
    lines = ['digraph ' + str(rsys.name) + '{']
    ind = '  '  # indentation

    def add_vertex(sn, num, reac):
        snum = str(num) if num > 1 else ''
        name = getattr(substances[sn], 'tex_name' if tex else 'name')
        lines.append(ind + '"{}" -> "{}" [label ="{}"];'.format(
            *((name, rid, snum) if reac else (rid, name, snum))
        ))

    for ri, rxn in enumerate(rsys.rxns):
        rid = rprefix + str(ri+rref0)
        lines.append(ind + '{')
        lines.append(ind*2 + 'node ' + nodeparams.format(rid))
        lines.append(ind*2 + rid)
        lines.append(ind + '}')
        for sn, num in rxn.reactants.items():
            add_vertex(sn, num, True)
        for sn, num in rxn.products.items():
            add_vertex(sn, num, False)
    lines.append('}')
    return lines


def rsys2graph(rsys, substances, outpath, prog=None, **kwargs):
    """
    Convenience function to call `rsys2dot` and write output to file
    and render the graph

    Parameters
    ==========
    rsys: ReactionSystem
    substances: sequence of Substance instances
    outpath: path to graph to be rendered
    prog: command to render DOT file (default: dot)
    **kwargs: parameters to pass along to `rsys2dot`

    Exapmles
    ========

    >>> rsys2graph(rsys, sbstncs, '/tmp/out.png')  # doctest: +SKIP
    """
    lines = rsys2dot(rsys, substances, **kwargs)
    if outpath.endswith('.tex'):
        cmds = [prog or 'dot2tex']
    else:
        cmds = [prog or 'dot', '-T'+outpath.split('.')[-1]]
    tmpfh = tempfile.NamedTemporaryFile('wt')
    tmpfh.writelines(lines)
    p = subprocess.Popen(cmds + [tmpfh.name, '-o', outpath])
    retcode = p.wait()
    if retcode:
        fmtstr = "{}\n returned with exit status {}"
        raise RuntimeError(fmtstr.format(' '.join(cmds), retcode))
