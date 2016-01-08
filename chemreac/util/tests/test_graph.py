# -*- coding: utf-8 -*-

import os
import shutil
import tempfile
from chemreac.chemistry import Reaction, ReactionSystem, mk_sn_dict_from_names
from chemreac.util.graph import rsys2dot, rsys2graph
from chemreac.util.testing import slow


def _get_rsys():
    sbstncs = mk_sn_dict_from_names('AB', latex_name=(r'$\mathcal{A}$',
                                                      r'$\mathcal{B}$'))
    r1 = Reaction({'A': 2}, {'B': 1}, k=3.0)
    rsys = ReactionSystem([r1], sbstncs)
    return rsys


def test_rsys2dot():
    rsys = _get_rsys()
    assert list(map(str.strip, rsys2dot(rsys))) == [
        'digraph None{',
        '{',
        'node [label=r1 shape=diamond]',
        'r1',
        '}',
        '"A" -> "r1" [label ="2"];',
        '"r1" -> "B" [label =""];',
        '}'
    ]


@slow
def test_rsys2graph():
    rsys, sbstncs = _get_rsys()
    tempdir = tempfile.mkdtemp()
    try:
        rsys2graph(rsys, sbstncs, os.path.join(tempdir, 'out.png'))
        rsys2graph(rsys, sbstncs, os.path.join(tempdir, 'out.ps'))
        rsys2graph(rsys, sbstncs, os.path.join(tempdir, 'out.tex'))
    finally:
        shutil.rmtree(tempdir)
