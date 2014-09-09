========
chemreac
========

.. image:: https://travis-ci.org/bjodah/chemreac.png?branch=master
   :target: https://travis-ci.org/bjodah/chemreac
   :alt: Build status
.. image:: https://readthedocs.org/projects/chemreac/badge/?version=latest
   :target: http://chemreac.readthedocs.org/
   :alt: Documentation Status
.. image:: https://coveralls.io/repos/bjodah/chemreac/badge.png?branch=master
   :target: https://coveralls.io/r/bjodah/chemreac?branch=master
   :alt: Test coverage

chremreac aims to be a library collecting tools and utilities for
modeling of chemical kinetics problems, primarily in the context of
aqueous phase with external radiation fields. 

Documentation
=============

the documentation is found at http://chemreac.readthedocs.org/

Installation
============
.. install-start

Below you will find instructions for installation. You may also
look in ``scripts/`` folder for automated install scripts used
in continuous integration.

Prerequisites
-------------

- C++ compiler with C++11 support (e.g. GCC >= 4.8)
- Fortran compiler with ISO_C_BINDING support (Fortran 2003 standard) (e.g. gfortran)
- LAPACK
- Sundials 2.5
- Python (2.7)
    
In addition to python, the following python packages are requiered:

- argh>=0.25.0
- numpy>=1.8.1
- cython>=0.19.1
- mako>=0.5.0
- quantities>=0.10.1
- pytest>=2.5.2
- scipy>=0.9.0
- matplotlib>=1.4.0
- periodictable>=1.4.1
- future>=0.12.3
- https://github.com/bjodah/pycompilation/archive/v0.3.3.tar.gz
- https://github.com/bjodah/pycodeexport/archive/v0.0.4.tar.gz
- https://github.com/sympy/sympy/archive/master.zip

For rendering the documentation you also need:

- `Sphinx <http://sphinx-doc.org/>`_
- `numpydoc <https://pypi.python.org/pypi/numpydoc>`_
- `sphinx_rtd_theme <https://pypi.python.org/pypi/sphinx_rtd_theme>`_

to run all the tests you also need:

- valgrind
- graphviz
- dot2tex

Once non-python prerequisites are in installed, you may procede e.g. as:

::

    $ git clone https://bitbucket.org/bjodah/chemreac.git
    $ cd chemreac
    $ pip install --user --upgrade -r requirements.txt
    $ python setup.py build_ext --inplace
    $ py.test


the above procedure works on Ubuntu 14.04 for example. 

To specify an alternative LAPACK lib set the environment variable LLAPACK, e.g.:

::

    $ LLAPACK=openblas python setup.py build_ext --inplace


Tests
-----
Run py.test, possibly with explicit PYTHONPATH (if build_ext --inplace was used)

::

    $ PYTHONPATH=`pwd`:$PYTHONPATH py.test

All tests should pass (or xfail).

.. install-end

Status
======

Continuous integration
----------------------
Build status at drone.io (uses OpenMP and tests sundials backend):

.. image:: https://drone.io/github.com/bjodah/chemreac/status.png
   :target: https://drone.io/github.com/bjodah/chemreac/latest

Build status at travis-ci (neither OpenMP nor sundials backend):

.. image:: https://travis-ci.org/bjodah/chemreac.png?branch=master
   :target: https://travis-ci.org/bjodah/chemreac


TODO
----
- Better defined "isolating" boundary conditions for:
    - logx
    - cylindrical/spherical
    - drift
    - and combinations of above.
- Look into logx+lrefl

License
=======
Open Source. Released under the very permissive "simplified
(2-clause) BSD license". See ``LICENSE.txt`` for further details.

Author
======
Björn Dahlgren, contact:
 - gmail adress: bjodah
 - kth.se adress: bda
