.. chemreac documentation master file, created by
   sphinx-quickstart on Tue Aug 26 17:03:58 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to chemreac's documentation!
====================================

chemreac is an open-source library for modeling chemical kinetcs in either:

- Homogeneous bulk solution
    - Kinetics governed by law of mass action
- Non-homogeneous solution where concentration profile has either:
    - planar symmetry
    - cylindrical symmetry
    - spherical symmetry

For the non-homogeneous case the user may choose:

- reflective or interpolating boundary conditions
- number of stencil points (3, 5 or 7)
- arbitrarily spaced grid
- calculate the electric field from concentrations for advection (drift).

Interfaces are provided to Sundials (CVODE) at the C++ level and to both 
Sundials and ODEPACK at the python level.

Currently the code is written with the following assumptions:

- isothermal conditions

TODO:

- Screening (Debye-Hückel)

Advection/Diffusion/Reaction models.
The model is formulated as the Smoluchovski equation and is discretied
in one dimension. The model is implmemnted in a C++ class with Python bindings
for ease of access.

Contents:

.. toctree::
   :maxdepth: 2

   core.rst
   integrate.rst
   chemistry.rst
   util.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

