# -*- coding: utf-8 -*-

"""
This file is an abstraction layer over ``quantities`` (makes it possible to
use another underlying package in the future).
"""

import quantities as pq

# Base units
metre = pq.metre
kilogram = pq.kilogram
ampere = pq.ampere
second = pq.second
Kelvin = pq.Kelvin
candela = pq.candela
mole = pq.mole

SI_base = {
    'length': metre,
    'mass': kilogram,
    'time': second,
    'current': ampere,
    'temperature': Kelvin,
    'luminous_intensity': candela,
    'amount': mole
}

# Convenience
joule = pq.joule
Gray = joule/kilogram
MeV = pq.MeV
centimetre = pq.centimetre
gram = pq.gram
dm = pq.UnitQuantity('decimeter',  pq.m / 10.0,  symbol='dm')
molar = pq.UnitQuantity('molar',  pq.mole / dm ** 3,  symbol='M')
perMolar_perSecond = 1/molar/pq.s
per100eV = pq.UnitQuantity('per_100_eV',  1/(100*pq.eV),  symbol='(100eV)**-1')
umol = pq.UnitQuantity('micromole',  pq.mole/1e6,  u_symbol=u'μmol')
umol_per_J = umol / pq.joule

# Constants
N_A = pq.constants.Avogadro_constant.rescale(1/pq.mol)
unitless_old_G_value_to_SI_G_value = (per100eV/N_A).rescale(umol_per_J)


# Utilities
def isunitless(expr):
    if hasattr(expr, 'dimensionality'):
        return expr.dimensionality == pq.dimensionless.dimensionality
    return True


def unitof(expr):
    try:
        return expr.units
    except AttributeError:
        return 1


def unit_registry_to_human_readable(unit_registry):
    new_registry = {}
    for k in SI_base:
        if unit_registry[k] is 1:
            new_registry[k] = 1, 1
        else:
            dim_list = list(unit_registry[k].dimensionality)
            if len(dim_list) != 1:
                raise TypeError("Compound units not allowed: {}".format(
                    dim_list))
            u_symbol = dim_list[0].u_symbol
            # u_symbol = unit_registry[k].u_symbol
            new_registry[k] = float(unit_registry[k]), u_symbol
    return new_registry


def unit_registry_from_human_readable(unit_registry):
    new_registry = {}
    for k in SI_base:
        factor, u_symbol = unit_registry[k]
        if u_symbol == 1:
            unit_quants = [1]
        else:
            unit_quants = list(pq.Quantity(0, u_symbol).dimensionality.keys())

        if len(unit_quants) != 1:
            raise TypeError("Unkown UnitQuantity: {}".format(unit_registry[k]))
        else:
            new_registry[k] = factor*unit_quants[0]
    return new_registry
