# -*- coding: utf-8 -*-

from ..constants import Boltzmann_constant, elementary_charge


def electrical_mobility_from_D(Dcoeff, charge, Temp, kB=None, elem_chg=None):
    """
    Calculates the electrical mobility through Einstein-Smoluchowski relation.

    Paraneters
    ----------
    Dcoeff: float
        Diffusion coefficient
    charge: integer
        charge of the species
    Temp: float
        Absolute temperature
    kB: float (optional)
        Boltzmann's constant. Defaults to value with units (which requires the
        other parameters to have consistent units)
    elem_chg: float (optional)
        elementary charge. Defaults to value with units (which requires the
        other parameters to have consistent units)

    Returns
    -------
    Electrical mobility (with units unless kB/elem_chg are overridden)
    """
    kB = Boltzmann_constant if kB is None else kB
    elem_chg = elementary_charge if elem_chg is None else elem_chg
    return Dcoeff*charge*elem_chg/(kB*Temp)
