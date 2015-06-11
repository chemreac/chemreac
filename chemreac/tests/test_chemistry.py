# -*- coding: utf-8 -*-

import quantities as pq
from periodictable import formula

from chemreac.chemistry import (
    molar, Substance, Reaction, ReactionSystem,
    mk_sn_dict_from_names
)


def test_to_ReactionDiffusion():
    sbstncs = mk_sn_dict_from_names('AB')
    r1 = Reaction({'A': 2}, {'B': 1}, k=3.0)
    rsys = ReactionSystem([r1])
    rd = rsys.to_ReactionDiffusion(sbstncs)
    assert rd.stoich_active == [[0, 0]]
    assert rd.stoich_prod == [[1]]
    assert rd.k == [3.0]


def test_Substance():
    formula_H2O = formula('H2O')
    H2O = Substance(name='H2O',  charge=0, formula=formula_H2O,
                    tex_name=r'$\mathrm{H_{2}O}$', pKa=14)
    OH_m = Substance(name='OH-',  charge=-1, formula=formula('OH'),
                     tex_name=r'$\mathrm{OH^{-}}$')
    assert sorted([OH_m, H2O]) == [H2O, OH_m]
