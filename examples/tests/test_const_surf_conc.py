# -*- coding: utf-8 -*-

from __future__ import (
    print_function, division, absolute_import, unicode_literals
)

from const_surf_conc import integrate_rd


def DISABLED_test_const_surf_conc():
    # It should be possible to get this to work again
    # but it will require some tweaking.

    tout, Cout, info, rd, tot_ave_rmsd = integrate_rd(
        D=2e-3, t0=1, tend=13, x0=1e-6, xend=1, N=1024, nt=42,
        logt=False, logy=True, logx=True, k=1.0, nstencil=3, scaling=1e-20,
        factor=1e12, atol=1e-8, rtol=1e-8)
    assert tot_ave_rmsd < 500
