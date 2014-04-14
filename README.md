chemreac
========
Chemical Reaction-Diffusion ODE systems discretized in one
dimension. Currently flat, spherical and cylindrical geometries are
implemented. Kinetics are governed by the law of mass action (with the
exception that it is allowed to make a difference between kinetically
active and stoichiometric reactants - e.g. solvent)

The Jacobian matrix can be evaluated and stored in several
formats. (currently dense row-major, dense col-major and banded
formats are supported).

The system of ODE is evaluated in the C++ class chemreac::ReactionDiffusion
(see `src/chemreac_template.cpp`). It is conveniently accessible from Python
(wrapped using Cython). For ease of use, a simple wrapper to
scipy.integrate.ode is provided in `chemreac.integrate.run`.

Setup
=====
```
    git clone https://bitbucket.org/bjodah/chemreac.git
    cd chemreac
    python pip install --upgrade -r requirements.txt
    python setup.py build_ext --inplace
```

should be enough (tested on Ubuntu 12.04). Add to $PYTHONPATH to use
or ``python setup.py install --user``

Tests
=====
Run py.test
``py.test``
(requires make, python-pytest), eventually source distributions will
include ``runtests.py`` which enables ``python setup.py test`` to work
without having pytest installed.

Prerequisites
=============
In addition to Python packages listed in ``requirements.txt`` you also need:
C++11 compliant compiler

TODO
======
- Jacobian generated for diffusion incorrect (with_jacobian=False performs MUCH better).
    * Need to derive test cases by hand.


License
=======
Open Source. Released under the very permissive "simplified
(2-clause) BSD license". See ``LICENSE.txt`` for further details.

Author
======
Björn Dahlgren, contact (gmail adress): bjodah
python analytic_diffusion.py -g s -N 3 --t0 10 --tend 100
python analytic_diffusion.py -g s --t0 10 --tend 100
