# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function)

from four_species import integrate_rd


def test_four_species():
    integr = integrate_rd(integrator='cvode', nt=2, dense_output=None)
    assert integr.info['success']
    assert len(integr.tout) > 2
