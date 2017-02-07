# -*- coding: utf-8 -*-
"""
Python package for modeling chemical kinetics with diffusion and drift.

chemreac is a python library for solving chemical kinetics problems
with possible diffusion and drift contributions. It is implemented
by solving the deterministic Smoluchovski equation for discretized
1D symmetric systems. Wrappers around SciPy's ODE integrators and
Sundials CVode package are provided.
"""
from __future__ import absolute_import, division, print_function

from ._release import __version__

from .core import ReactionDiffusion, Geom_names


def get_include():
    from pkg_resources import resource_filename, Requirement
    return resource_filename(Requirement.parse(__name__),
                             '%s/include' % __name__)
