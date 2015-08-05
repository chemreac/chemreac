from chemreac.util.physchem import electrical_mobility_from_D
from chemreac.units import metre, second, kelvin, coulomb, kilogram, allclose


def test_electrical_mobility_from_D():
    D = 3*metre/second**2
    z = -2
    T = 100*kelvin
    mu = electrical_mobility_from_D(D, z, T)
    e = 1.60217657e-19 * coulomb
    kB = 1.3806488e-23 * metre**2 * kilogram * second**-2 / kelvin
    ref = z*e/kB/T*D
    assert allclose(mu, ref, rtol=1e-5)
