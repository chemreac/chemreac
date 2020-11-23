import numpy as np

from steady_state_approx import integrate_rd


def _assert_mp(t, ydot, yout):
    from mpmath import odefun
    f = odefun(ydot, 0, [1, 0, 0])
    for tidx, tval in enumerate(t):
        mpvals = f(tval)
        for sidx in range(3):
            assert abs(mpvals[sidx] - yout[tidx, 0, sidx]) < 1e-8


def test_Bss_approxiamtion():
    t, yout, A_ssB, A_ssB_2fast, ydot = integrate_rd(
        1.0, 1e-3, 220, 217)
    assert np.allclose(yout[:, 0, 0], A_ssB)
    _assert_mp(t, ydot, yout)


def test_Bss_2fast_approxiamtion():
    t, yout, A_ssB, A_ssB_2fast, ydot = integrate_rd(
        1.0, 1e-10, 220e5, 217)
    assert np.allclose(yout[:, 0, 0], A_ssB_2fast)
    # Too long runtime:
    # _assert_mp(t, ydot, yout)
