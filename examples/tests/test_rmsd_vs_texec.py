# -*- coding: utf-8 -*-

from collections import OrderedDict

from rmsd_vs_texec import main


def test_integrate_rd():
    main(OrderedDict([
        ('nstencil', [3, 5]),
        ('N', [40, 80]),
        ('method', ['bdf', 'adams'])
    ]))
