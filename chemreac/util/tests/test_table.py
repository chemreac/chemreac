# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

import shutil
import tempfile

import pytest
from chempy import Reaction, Substance, ReactionSystem
from chemreac import ReactionDiffusion
from chemreac.util.table import radyields2pdf_table


def _get_rsys():
    r1 = Reaction({'A': 2}, {'B': 1}, 3.0)
    A = Substance('A', latex_name='\\ensuremath{\\boldsymbol{A}}')
    B = Substance('B', latex_name='\\ensuremath{\\boldsymbol{B}}')
    rsys = ReactionSystem([r1], [A, B])
    return rsys


@pytest.mark.xfail
def test_radyields2pdf_table():
    rsys = _get_rsys()
    rd = ReactionDiffusion.from_ReactionSystem(rsys)
    tempdir = tempfile.mkdtemp()
    try:
        radyields2pdf_table(rd, tempdir)
    finally:
        shutil.rmtree(tempdir)
