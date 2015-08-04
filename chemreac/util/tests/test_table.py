# -*- coding: utf-8 -*-

import shutil
import tempfile

import pytest


from chemreac.util.table import (
    rsys2tablines, rsys2table, rsys2pdf_table, radyields2pdf_table
)
from .test_graph import _get_rsys


def test_rsys2tablines():
    assert rsys2tablines(*_get_rsys(), tex=False) == [
        '1 & 2A & $\\rightarrow$ & B & 3 & 1 & None'
    ]


def test_rsys2table():
    assert rsys2table(*_get_rsys()) == r"""
\begin{table}
\centering
\label{tab:none}
\caption[None]{None}
\begin{tabular}{lllllll}
\toprule
Id. & Reactants &  & Products & {Rate constant} & Unit & Ref \\
\midrule
1 & 2$\mathcal{A}$ & $\rightarrow$ & $\mathcal{B}$ & 3 & 1 & None \\
\bottomrule
\end{tabular}
\end{table}"""


@pytest.mark.parametrize('longtable', (True, False))
def test_rsys2pdf_table(longtable):
    rsys, sbstncs = _get_rsys()
    tempdir = tempfile.mkdtemp()
    try:
        rsys2pdf_table(rsys, sbstncs, tempdir, longtable=longtable)
    finally:
        shutil.rmtree(tempdir)


def test_rsys2pdf_table_no_output_dir():
    rsys, sbstncs = _get_rsys()
    rsys2pdf_table(rsys, sbstncs, save=False)


def test_radyields2pdf_table():
    rsys, sbstncs = _get_rsys()
    rd = rsys.to_ReactionDiffusion(sbstncs)
    tempdir = tempfile.mkdtemp()
    try:
        radyields2pdf_table(rd, tempdir)
    finally:
        shutil.rmtree(tempdir)
