# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

from itertools import product
import sys

import pytest

from chemreac.units import molar, allclose, SI_base_registry, umol, decimetre

import with_units

TR_FLS = (True, False)


@pytest.mark.skipif(sys.version_info[0] > 2, reason='pytest/python3 incompatible with this test')
@pytest.mark.parametrize('params', list(product(TR_FLS, TR_FLS)))
def test_simple(params):
    logy, logt = params
    integr, Cref, rd = with_units.main(logy=logy, logt=logt)
    Cout = integr.with_units('Cout')
    if logy or logt:
        assert allclose(Cout[:, 0, :], Cref, atol=2e-5*molar)
    else:
        assert allclose(Cout[:, 0, :], Cref, atol=1e-6*molar)


@pytest.mark.skipif(sys.version_info[0] > 2, reason='pytest/python3 incompatible with this test')
def test_simple_other_units():
    unit_registry = SI_base_registry.copy()
    unit_registry['length'] = decimetre
    unit_registry['amount'] = umol
    integr, Cref, rd = with_units.main(unit_registry=unit_registry)
    Cout = integr.with_units('Cout')
    assert allclose(Cout[:, 0, :], Cref, atol=1e-6*molar)
