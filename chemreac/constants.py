# -*- coding: utf-8 -*-

import quantities as pq

from .units import per100eV, umol_per_J

faraday = pq.constants.Faraday_constant
vacuum_permittivity = pq.constants.vacuum_permittivity
Avogadro_constant = pq.constants.Avogadro_constant.rescale(1/pq.mol)

# Conversion factors
unitless_old_G_value_to_SI_G_value = (per100eV/Avogadro_constant).rescale(
    umol_per_J)
