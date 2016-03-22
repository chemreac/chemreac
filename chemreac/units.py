# -*- coding: utf-8 -*-

"""
This file is an abstraction layer over ``quantities`` (makes it possible to
use another underlying package in the future).
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import quantities as pq

from chempy.units import (
    allclose, default_units, get_derived_unit, is_unitless, linspace,
    SI_base_registry, to_unitless, unit_of, unit_registry_to_human_readable,
    unit_registry_from_human_readable
)


metre = default_units.metre
kilogram = default_units.kilogram
ampere = default_units.ampere
second = default_units.second
kelvin = default_units.kelvin
candela = default_units.candela
mole = default_units.mole
coulomb = default_units.coulomb
joule = default_units.joule
gray = default_units.gray
eV = default_units.eV
MeV = default_units.MeV
metre = default_units.metre
decimetre = default_units.decimetre
centimetre = default_units.centimetre
micrometre = default_units.micrometre
nanometre = default_units.nanometre
gram = default_units.gram
molar = default_units.molar
hour = default_units.hour
day = default_units.day
perMolar_perSecond = default_units.perMolar_perSecond
per100eV = default_units.per100eV
umol = default_units.umol
umol_per_J = default_units.umol_per_J
