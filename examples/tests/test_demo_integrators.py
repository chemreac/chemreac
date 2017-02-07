# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import demo_integrators as di


def test_demo_integrators_ff():
    di.main(False, False, 1)


def test_demo_integrators_ft():
    di.main(False, True, 10000, 1034*1.05)


def test_demo_integrators_tf():
    di.main(True, False, 8)


def test_demo_integrators_tt():
    di.main(True, True, 9000, 1151*1.05)
