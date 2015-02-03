# -*- coding: utf-8 -*-

import shutil
import tempfile

import pytest


from chemreac.util.table import rsys2tablines, rsys2table, rsys2pdf_table
from .test_graph import _get_rsys


def test_rsys2tablines():
    assert rsys2tablines(*_get_rsys(), tex=False) == [
        '1 & 2A & $\\rightarrow$ & B & 3 & None'
    ]


def test_rsys2table():
    assert rsys2table(*_get_rsys()) == r"""
\begin{table}
\centering
\label{tab:none}
\caption[None]{None}
\begin{tabular}{llllll}
\toprule
Id. & Reactants &  & Products & Rate constant & Ref \\
\midrule
1 & 2$\mathcal{A}$ & $\rightarrow$ & $\mathcal{B}$ & 3 & None \\
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
