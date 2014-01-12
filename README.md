chemreac
========
Chemical Reaction-Diffusion ODE systems discretized in one
dimension. Currently flat and spherical geometry is
implemented. Kinetics are governed by the law of mass action (with the
exception that we allow to make a difference between kinetically
active and stoichiometric reactants - e.g. solvent may participate)

The Jacobian can be evaluated and stored in several
formats. (currently dense row-major, dense col-major and banded
formats are supported).

The ODE systems is computed in the C++ class ReactionDiffusion
(`cpp_chem_template.cpp`). It is conveniently accessible from Python
(wrapped using Cython). A simple wrapper to scipy.integrate.ode is
provided in `chemreac.integrate.run`

Setup
=====
``` python setup.py build_ext --inplace ``` is enough on most *NIX machines.

Tests
=====
Run `./py.test`

Prerequisites
=============
C++11 compliant compiler
Python
Cython
Pycompilation 0.2.4

License
=======
Open Soucrce. Released under the very permissive "simplified
(2-clause) BSD license". See LICENCE.txt for further details.

Author
======
Björn Dahlgren, contact (gmail adress): bjodah
python -m pudb /home/bjorn/.local/bin/xdress --debug
env DISTUTILS_DEBUG=1 python -m pudb setup.py build
