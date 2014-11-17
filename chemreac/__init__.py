# -*- coding: utf-8 -*-
"""
Python extension for reaction diffusion.
"""

import os
import subprocess

from .release import __version__

from .core import (
    ReactionDiffusion, DENSE, BANDED, SPARSE, FLAT,
    CYLINDRICAL, SPHERICAL, Geom_names
)
