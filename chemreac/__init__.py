# -*- coding: utf-8 -*-
"""
Python extension for reaction diffusion.
"""

from .release import __version__

from .core import (
    ReactionDiffusion, DENSE, BANDED, SPARSE, FLAT,
    CYLINDRICAL, SPHERICAL, Geom_names
)
