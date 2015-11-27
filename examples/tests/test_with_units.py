# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from future.builtins import *

from itertools import product

import pytest

from chemreac.units import molar, allclose, SI_base, umol, decimetre

import with_units

TR_FLS = (True, False)


@pytest.mark.parametrize('params', list(product(TR_FLS, TR_FLS)))
def test_simple(params):
    logy, logt = params
    integr, Cref, rd = with_units.main(logy=logy, logt=logt)
    if logy or logt:
        assert allclose(integr.Cout[:, 0, :], Cref, atol=2e-5*molar)
    else:
        assert allclose(integr.Cout[:, 0, :], Cref, atol=1e-6*molar)


def test_simple_other_units():
    unit_registry = SI_base.copy()
    unit_registry['length'] = decimetre
    unit_registry['amount'] = umol
    integr, Cref, rd = with_units.main(unit_registry=unit_registry)
    assert allclose(integr.Cout[:, 0, :], Cref, atol=1e-6*molar)
