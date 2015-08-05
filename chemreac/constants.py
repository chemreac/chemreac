# -*- coding: utf-8 -*-

import quantities as pq

from .units import to_unitless, get_derived_unit

Faraday_constant = pq.constants.Faraday_constant
vacuum_permittivity = pq.constants.vacuum_permittivity
Avogadro_constant = pq.constants.Avogadro_constant.rescale(1/pq.mol)
Boltzmann_constant = pq.constants.Boltzmann_constant
elementary_charge = pq.constants.elementary_charge


def get_unitless_constant(registry, name):
    named_const = {
        'Faraday_constant': Faraday_constant,
        'vacuum_permittivity': vacuum_permittivity,
    }
    if registry is None:
        return float(named_const[name].definition)
    named_const_units = {
        'Faraday_constant': (get_derived_unit(registry, 'charge') /
                             get_derived_unit(registry, 'amount')),
        'vacuum_permittivity': get_derived_unit(registry, 'permittivity')
    }
    return to_unitless(named_const[name], named_const_units[name])
