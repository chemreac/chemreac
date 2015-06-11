# -*- coding: utf-8 -*-

"""
This file is an abstraction layer over ``quantities`` (makes it possible to
use another underlying package in the future).
"""
from __future__ import absolute_import, division, print_function

import numpy as np
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


def get_derived_unit(registry, name):
    if registry is None:
        return 1.0
    derived = {
        'diffusion': registry['length']**2/registry['time'],
        'electrical_mobility': (registry['current']*registry['time']**2 /
                                registry['mass']),
        'permittivity': (registry['current']**2*registry['time']**4 /
                         (registry['length']**3*registry['mass'])),
        'charge': registry['current']*registry['time'],
        'energy': registry['mass']*registry['length']**2/registry['time']**2,
        'concentration': registry['amount']/registry['length']**3,
        'density': registry['mass']/registry['length']**3,
    }
    derived['radiolytic_yield'] = registry['amount']/derived['energy']
    try:
        return derived[name]
    except KeyError:
        return registry[name]


# Convenience
joule = pq.joule
gray = pq.gray
eV = pq.eV
MeV = pq.MeV
metre = pq.metre
decimetre = dm = pq.UnitQuantity('decimetre',  pq.m / 10.0,  u_symbol='dm')
centimetre = pq.centimetre
micrometre = pq.micrometre
nanometre = pq.nanometre
gram = pq.gram
molar = pq.UnitQuantity('molar',  pq.mole / dm ** 3,  u_symbol='M')
hour = pq.hour
perMolar_perSecond = 1/molar/pq.s
per100eV = pq.UnitQuantity('per_100_eV',
                           1/(100*pq.eV*pq.constants.Avogadro_constant),
                           u_symbol='(100eV)**-1')
umol = pq.UnitQuantity('micromole',  pq.mole/1e6,  u_symbol=u'μmol')
umol_per_J = umol / pq.joule


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


def to_unitless(value, new_unit=None):
    if new_unit is None:
        new_unit = pq.dimensionless
    if isinstance(value, (list, tuple)):
        return np.array([to_unitless(elem, new_unit) for elem in value])
    try:
        result = (value*pq.dimensionless/new_unit).rescale(pq.dimensionless)
        if result.ndim == 0:
            return float(result)
        else:
            return np.asarray(result)
    except TypeError:
        return np.array([to_unitless(elem, new_unit) for elem in value])


def rescale(value, new_unit):
    if isinstance(new_unit, pq.quantity.Quantity):
        return (value*pq.dimensionless/new_unit).rescale(
            pq.dimensionless)*new_unit
    else:
        if isinstance(value, pq.quantity.Quantity):
            raise ValueError("Cannot rescale {} to dimensionless".format(
                value))
        else:
            try:
                return value/new_unit
            except TypeError:
                return [elem/new_unit for elem in value]


def unit_registry_to_human_readable(unit_registry):
    if unit_registry is None:
        return None
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
    if unit_registry is None:
        return None
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


def allclose(a, b, rtol=1e-8, atol=None):
    d = np.abs(a - b)
    lim = np.abs(a)*rtol
    if atol is not None:
        lim += atol
    return np.all(d < lim)
