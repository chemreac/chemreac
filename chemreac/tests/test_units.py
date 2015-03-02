from collections import defaultdict

from chemreac.units import (
    isunitless, unitof, unit_registry_to_human_readable,
    unit_registry_from_human_readable, SI_base,
    metre, kilogram, second, ampere, Kelvin, candela, mole,
    dm
)


def test_isunitless():
    assert not isunitless(1*dm)
    assert isunitless(1)


def test_unitof():
    assert unitof(0.1*mole/dm**3) == mole/dm**3
    assert unitof(7) == 1


def test_unit_registry_to_human_readable():
    # Not as much human readable as JSON serializable...
    d = defaultdict(lambda: 1)
    assert unit_registry_to_human_readable(d) == dict(
        (x, (1, 1)) for x in SI_base.keys())

    ur = {
        'length': 1e3*metre,
        'mass': 1e-2*kilogram,
        'time': 1e4*second,
        'current': 1e-1*ampere,
        'temperature': 1e1*Kelvin,
        'luminous_intensity': 1e-3*candela,
        'amount': 1e4*mole
    }
    assert unit_registry_to_human_readable(ur) == {
        'length': (1e3, 'm'),
        'mass': (1e-2, 'kg'),
        'time': (1e4, 's'),
        'current': (1e-1, 'A'),
        'temperature': (1e1, 'K'),
        'luminous_intensity': (1e-3, 'cd'),
        'amount': (1e4, 'mol')
    }
    assert unit_registry_to_human_readable(ur) != {
        'length': (1e2, 'm'),
        'mass': (1e-2, 'kg'),
        'time': (1e4, 's'),
        'current': (1e-1, 'A'),
        'temperature': (1e1, 'K'),
        'luminous_intensity': (1e-3, 'cd'),
        'amount': (1e4, 'mol')
    }


def test_unit_registry_from_human_readable():
    hr = unit_registry_to_human_readable(defaultdict(lambda: 1))
    assert hr == dict((x, (1, 1)) for x in SI_base.keys())
    ur = unit_registry_from_human_readable(hr)
    assert ur == dict((x, 1) for x in SI_base.keys())

    hr = unit_registry_to_human_readable(SI_base)
    assert hr == {
        'length': (1.0, 'm'),
        'mass': (1.0, 'kg'),
        'time': (1.0, 's'),
        'current': (1.0, 'A'),
        'temperature': (1.0, 'K'),
        'luminous_intensity': (1.0, 'cd'),
        'amount': (1.0, 'mol')
    }
    ur = unit_registry_from_human_readable(hr)
    assert ur == SI_base

    ur = unit_registry_from_human_readable({
        'length': (1.0, 'm'),
        'mass': (1.0, 'kg'),
        'time': (1.0, 's'),
        'current': (1.0, 'A'),
        'temperature': (1.0, 'K'),
        'luminous_intensity': (1.0, 'cd'),
        'amount': (1.0, 'mol')
    })
    assert ur == {
        'length': metre,
        'mass': kilogram,
        'time': second,
        'current': ampere,
        'temperature': Kelvin,
        'luminous_intensity': candela,
        'amount': mole
    }

    ur = unit_registry_from_human_readable({
        'length': (1e3, 'm'),
        'mass': (1e-2, 'kg'),
        'time': (1e4, 's'),
        'current': (1e-1, 'A'),
        'temperature': (1e1, 'K'),
        'luminous_intensity': (1e-3, 'cd'),
        'amount': (1e4, 'mol')
    })
    assert ur == {
        'length': 1e3*metre,
        'mass': 1e-2*kilogram,
        'time': 1e4*second,
        'current': 1e-1*ampere,
        'temperature': 1e1*Kelvin,
        'luminous_intensity': 1e-3*candela,
        'amount': 1e4*mole
    }

    assert ur != {
        'length': 1e2*metre,
        'mass': 1e-3*kilogram,
        'time': 1e2*second,
        'current': 1e-2*ampere,
        'temperature': 1e0*Kelvin,
        'luminous_intensity': 1e-2*candela,
        'amount': 1e3*mole
    }
