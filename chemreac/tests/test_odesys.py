# -*- coding: utf-8 -*-
from collections import defaultdict
import numpy as np
from chemreac import ReactionDiffusion
from chempy import ReactionSystem
from chempy.units import to_unitless, SI_base_registry, get_derived_unit, default_units as u

analytic = [
    lambda y0, k, t: (
        y0[0] * np.exp(-k[0]*t)),
    lambda y0, k, t: (
        y0[1] * np.exp(-k[1] * t) + y0[0] * k[0] / (k[1] - k[0]) *
        (np.exp(-k[0]*t) - np.exp(-k[1]*t))),
    lambda y0, k, t: (
        y0[2] + y0[1] * k[1] / (-k[1]) *
        (np.exp(-k[1]*t) - 1) +
        k[1] * k[0] * y0[0] / (k[1] - k[0]) *
        (1 / (-k[0]) * (np.exp(-k[0]*t) - 1) -
         1 / (-k[1]) * (np.exp(-k[1]*t) - 1)))
]


def _get_odesys(k):
    names = ['A', 'B', 'C']
    pns = ['kA', 'kB']
    rd = ReactionDiffusion(len(names), [[0], [1]], [[1], [2]], k=k,
                           substance_names=names, param_names=pns)
    return rd._as_odesys(k_from_params=lambda p: [p[k] for k in pns])


def test_decay():
    kA = 0.13
    odesys = _get_odesys([kA, 0])
    y0 = dict(A=3, B=1, C=0)
    t0, tend, nt = 5.0, 17.0, 42
    tout = np.linspace(t0, tend, nt+1)
    result = odesys.integrate(tout, y0)
    yref = np.array([y0['A']*np.exp(-kA*(tout-t0)),
                     y0['B']+y0['A']*(1-np.exp(-kA*(tout-t0)))]).transpose()
    assert np.allclose(result.yout[:, :2], yref)


def test_decay_params():
    odesys = _get_odesys([0, 0])
    y0 = 42, 7, 4
    k = .7, .3
    ic = dict(zip(odesys.names, y0))
    p = dict(zip('kA kB'.split(), k))
    tout, yout, info = odesys.integrate([0, 5], ic, p)
    yref = np.array([a(y0, k, tout) for a in analytic]).transpose()
    assert np.allclose(yout, yref)


def test_chained_parameter_variation():
    # A -> B
    names = ['A', 'B']
    rd = ReactionDiffusion(len(names), [], [], k=[],
                           substance_names=names, g_value_parents=[0], g_values=[[0, 1]])
    durations = [1., 3., 2.]
    y0 = [13., 7.]
    ic = dict(zip(names, y0))
    doserates = [.3, .11, .7]
    npoints = 3
    odesys = rd._as_odesys()
    res = odesys.chained_parameter_variation(
        durations, ic, {'doserate': doserates}, {'density': 1.0}, npoints=npoints)
    assert res.xout.size == npoints*len(durations) + 1
    assert res.xout[0] == 0
    assert np.all(res.yout[0, :] == y0)
    cumulative = 0.0
    for dr, dur in zip(doserates, durations):
        mask = (cumulative <= res.xout) & (res.xout <= cumulative + dur)
        cumulative += dur
        t, y = res.xout[mask], res.yout[mask, :]
        a, b = y[:, 0], y[:, 1]
        refa = y0[0]
        refb = y0[1] + t*dr*y0[0]
        assert np.allclose(refa, a)
        assert np.allclose(refb, b)


def test_chained_parameter_variation_from_ReactionSystem():
    g_E_mol_J = 2.1e-7
    rsys = ReactionSystem.from_string(
        """
        (H2O) -> e-(aq) + H+ + OH; Radiology(%.2e*mol/Gy)
        2 OH -> H2O2; 3.6e-9/M/s
        H+ + OH- -> H2O; 1.4e11/M/s
        H2O -> H+ + OH-; 1.4e-3/s
        N2O + e-(aq) -> N2 + O-; 9.6e9/M/s
        O- + H+ -> OH; 1e11/M/s
        """ % g_E_mol_J  # neglecting a large body of reactions (just a test-case after all)
    )
    ureg = SI_base_registry
    field_u = get_derived_unit(ureg, 'doserate') * get_derived_unit(ureg, 'density')
    rd = ReactionDiffusion.from_ReactionSystem(rsys, fields=[0*field_u], unit_registry=ureg)
    odesys = rd._as_odesys()
    npoints = 5
    durations = [59*u.second, 42*u.minute, 2*u.hour]
    doserates = [135*u.Gy/u.s, 11*u.Gy/u.s, 180*u.Gy/u.minute]
    M = u.molar
    ic = defaultdict(lambda: 0*M, {'H2O': 55.4*M, 'H+': 1e-7*M, 'OH-': 1e-7*M, 'N2O': 20e-3*M})
    dens_kg_dm3 = 0.998

    result = odesys.chained_parameter_variation(
        durations, ic, {'doserate': doserates}, {'density': dens_kg_dm3*1e3*u.kg/u.m**3},
        npoints=npoints)

    N2_M = to_unitless(result.named_dep('N2'), u.M)
    H2O2_M = to_unitless(result.named_dep('H2O2'), u.M)

    e_accum_molar = 0
    for i, (dur, dr) in enumerate(zip(durations, doserates)):
        dur_s = to_unitless(dur, u.s)
        dr_Gy_s = to_unitless(dr, u.Gy/u.s)
        local_ts = np.linspace(0, dur_s, npoints+1)
        # local_ic = {k: result.named_dep(k)[i*npoints] for k in odesys.names}
        for j, (lt, ld) in enumerate(zip(local_ts[1:], np.diff(local_ts))):
            e_accum_molar += ld*g_E_mol_J*dr_Gy_s*dens_kg_dm3
            assert abs(N2_M[i*npoints + j] - e_accum_molar)/e_accum_molar < 1e-3
            assert abs(H2O2_M[i*npoints + j] - 2*e_accum_molar)/e_accum_molar < 1e-3
