# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from future.builtins import *

from itertools import product

import pytest

from chemreac.units import molar, allclose, SI_base, umol, decimetre

import with_units

TR_FLS = (True, False)


@pytest.mark.parametrize('params', list(product(TR_FLS, TR_FLS)))
def test_simple(params):
    logy, logt = params
    integr, Cref = with_units.main(logy=logy, logt=logt)
    if logy or logt:
        assert allclose(integr.Cout[:, 0, :], Cref, atol=2e-5*molar)
    else:
        assert allclose(integr.Cout[:, 0, :], Cref, atol=1e-6*molar)


def test_simple_other_units():
    units = SI_base.copy()
    units['length'] = decimetre
    units['amount'] = umol
    integr, Cref = with_units.main(units=units)
    assert allclose(integr.Cout[:, 0, :], Cref, atol=1e-6*molar)
