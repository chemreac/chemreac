# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function

from itertools import product

import numpy as np
import pytest

from chemreac.util.testing import veryslow, slow
from auto_efield import integrate_rd


TR_FLS = (True, False)

COMBOS = list(product(
    TR_FLS, TR_FLS, [False], [False], [3]))

EXTRA_COMBOS = list(product(  # "py.test --veryslow" to invoke
    TR_FLS, TR_FLS, TR_FLS, TR_FLS, [3]))  # nstencil > 3 too slow..


def _test_perfect_overlap(params):
    (ly, lt, lr, rr, ns), geom, lx = params
    kwargs = {'geom': geom, 'logt': lt, 'logy': ly, 'logx': lx,
              'nstencil': ns, 'lrefl': lr, 'rrefl': rr}

    # zero offset (offset = 0.0)
    tout, Cout, info, rd = integrate_rd(D=1e-5, N=64, offset=0,
                                        sigma_q=13, **kwargs)
    assert info['success']
    delta_C = Cout[:, :, 0] - Cout[:, :, 1]

    # test that the species do not deviate (weak test)
    assert np.all(np.abs(delta_C) < 1e-3)


@pytest.mark.parametrize('params', list(product(COMBOS, 'fcs', TR_FLS)))
def test_perfect_overlap(params):
    _test_perfect_overlap(params)


@veryslow
@pytest.mark.parametrize('params', list(product(EXTRA_COMBOS, 'fcs', TR_FLS)))
def test_perfect_overlap_extended(params):
    _test_perfect_overlap(params)


# Mass convservation for a pair of gaussians at a offset
def _test_mass_conservation(params):
    (ly, lt, lr, rr, ns), geom, lx = params
    kwargs = {'geom': geom, 'logt': lt, 'logy': ly, 'logx': lx,
              'nstencil': ns, 'lrefl': lr, 'rrefl': rr}

    # separated gaussians (offset = 0.25)
    N = 64*(4 if ly else 1)*(8 if ns > 3 else 1)
    tout, Cout, info, rd = integrate_rd(
        D=-3e-8, t0=1e-6, tend=7, x0=0.1, xend=1.0, N=N,
        offset=.25, nt=25, sigma_q=101, **kwargs)
    assert info['success']
    mass_consv = np.array([
        [rd.integrated_conc(Cout[j, :, i]) for i in range(rd.n)]
        for j in range(tout.size)
    ])
    print(mass_consv - mass_consv[[0], :])
    assert np.all(np.abs(mass_consv - mass_consv[[0], :]) < 1e-5)
    # For current parameters gaussians are well separated at end
    # of simulation, hence we expect maximum Efield to have same
    # value as in the beginning of simulation
    efield_i = rd.calc_efield(Cout[0, :, :].flatten())
    efield_f = rd.calc_efield(Cout[-1, :, :].flatten())
    assert abs(np.max(np.abs(efield_i)) - np.max(np.abs(efield_f))) < 1e-6


@pytest.mark.parametrize('params', list(product(COMBOS, 'f', [False])))
def test_mass_conservation_flat(params):
    _test_mass_conservation(params)


@veryslow
@pytest.mark.parametrize('params', list(product(EXTRA_COMBOS, 'f', [False])))
def test_mass_conservation_flat_extended(params):
    _test_mass_conservation(params)


@pytest.mark.parametrize('params', list(product(COMBOS, 'cs', [False])))
# @pytest.mark.xfail
def test_mass_conservation_cyl_sph(params):
    _test_mass_conservation(params)


# Pair of half-gaussians centered at x0 (different sigma), using logx
def _test_pair_centered_at_x0_different_sigma(params):
    (ly, lt, lr, rr, ns), geom, lx = params
    kwargs = {'geom': geom, 'logt': lt, 'logy': ly, 'logx': lx,
              'nstencil': ns, 'lrefl': lr, 'rrefl': rr}

    # separated gaussians (offset = 0.25)
    N = 64*(8 if ly else 1)*(4 if ns > 3 else 1)
    tout, Cout, info, rd = integrate_rd(
        D=-3e-8, t0=1e-6, tend=7, x0=1e-6, xend=1.0, N=N,
        base=0, offset=0, nt=25, sigma_q=101,
        sigma_skew=0.1, **kwargs)
    assert info['success']
    mass_consv = np.array([
        [rd.integrated_conc(Cout[j, :, i]) for i in range(rd.n)]
        for j in range(tout.size)
    ])
    assert np.all(np.abs(mass_consv - mass_consv[[0], :]) < 1e-4)


@pytest.mark.parametrize('params', list(product(COMBOS, 'f', [True])))
def test_pair_centered_at_x0_different_sigma_flat_logx(params):
    _test_pair_centered_at_x0_different_sigma(params)


@slow
@pytest.mark.parametrize('params', list(product(EXTRA_COMBOS, 'f', [True])))
def test_pair_centered_at_x0_different_sigma_flat_logx_extended(params):
    _test_pair_centered_at_x0_different_sigma(params)


@slow
@pytest.mark.parametrize('params', list(product(COMBOS, 'cs', [True])))
@pytest.mark.xfail
def test_pair_centered_at_x0_different_sigma_cyl_sph_logx(params):
    _test_pair_centered_at_x0_different_sigma(params)
