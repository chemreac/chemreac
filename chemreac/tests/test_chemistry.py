# -*- coding: utf-8 -*-

import quantities as pq
from periodictable import formula

from chemreac.chemistry import molar, Substance, Henry


def test_Henry():
    h       = Henry(1.3e-3 * molar / pq.atm,
              1500 * pq.kelvin)
    T       = 293.15 * pq.kelvin
    ref_k_H = 0.0014164786921549344 * molar / pq.atm
    ref_c   = 0.0014164786921549344 * molar
    ref_P   = 101.325e3 * pq.pascal
    assert h.get_k_H_at_T(T) - ref_k_H < 1e-6 * molar / pq.atm
    assert h.get_c_at_T_and_P(T, ref_P) - ref_c < 1e-6 * molar
    assert h.get_P_at_T_and_c(T, ref_c) - ref_P < 1 * pq.pascal

def test_Substance():
    formula_H2O = formula('H2O')
    H2O = Substance(name = 'H2O',  charge = 0, formula = formula_H2O,
                  tex_name = r'$\mathrm{H_{2}O}$', pKa = 14)
    OH_m = Substance(name = 'OH-',  charge =-1, formula = formula('OH'),
                     tex_name = r'$\mathrm{OH^{-}}$')
    assert sorted([OH_m, H2O]) == [H2O, OH_m]
    assert str(H2O) == 'H2O'
    assert str(OH_m) == 'OH-'
    print(repr(H2O))
    assert repr(H2O) == "Substance('H2O', 0, {}, {}, '{}', 14, None, 0)".format(
        formula_H2O.mass, repr(formula_H2O), '$\mathrm{H_{2}O}$')
