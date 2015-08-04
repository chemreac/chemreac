# -*- coding: utf-8 -*-

import pytest
import quantities as pq
from periodictable import formula

from chemreac import ReactionDiffusion
from chemreac.units import molar, second
from chemreac.chemistry import (
    molar, Substance, Reaction, ReactionSystem,
    mk_sn_dict_from_names
)


def test_ReactionSystem__to_ReactionDiffusion():
    sbstncs = mk_sn_dict_from_names('AB')
    r1 = Reaction({'A': 2}, {'B': 1}, k=3.0)
    rsys = ReactionSystem([r1])
    rd = rsys.to_ReactionDiffusion(sbstncs)
    assert rd.stoich_active == [[0, 0]]
    assert rd.stoich_prod == [[1]]
    assert rd.k == [3.0]


def test_ReactionSystem__from_ReactionDiffusion():
    rd = ReactionDiffusion(2, [[0]], [[1]], [1])
    rsys = ReactionSystem.from_ReactionDiffusion(rd)
    assert len(rsys.rxns) == 1


def test_ReactionSystem():
    pass


def test_Substance():
    formula_H2O = formula('H2O')
    H2O = Substance(name='H2O',  charge=0, formula=formula_H2O,
                    tex_name=r'$\mathrm{H_{2}O}$', pKa=14)
    OH_m = Substance(name='OH-',  charge=-1, formula=formula('OH'),
                     tex_name=r'$\mathrm{OH^{-}}$')
    assert sorted([OH_m, H2O]) == [H2O, OH_m]


@pytest.mark.parametrize('equilibrium', (True, False))
def test_Reaction(equilibrium):
    rMrs = 1/molar/second
    formula_H2O = formula('H2O')
    H2O = Substance(name='H2O',  charge=0, formula=formula_H2O,
                    tex_name=r'$\mathrm{H_{2}O}$', pKa=14)
    H_p = Substance(name='H+',  charge=1, formula=formula('H'),
                    tex_name=r'$\mathrm{H^{+}}$')
    OH_m = Substance(name='OH-',  charge=-1, formula=formula('OH'),
                     tex_name=r'$\mathrm{OH^{-}}$')
    r1 = Reaction({H_p: 1, OH_m: 1}, {H2O: 1}, k=1.4e11*rMrs)
    r1_str = str(r1)
    for fragment in ('H+', 'H2O', 'OH-', '->'):
        assert fragment in r1_str
    r1_tex = r1.render(tex=True, equilibrium=equilibrium)
    for fragment in (H2O.tex_name, H_p.tex_name, OH_m.tex_name,
                     r'$\rightleftharpoons$' if equilibrium else
                     r'$\rightarrow$'):
        assert fragment in r1_tex
