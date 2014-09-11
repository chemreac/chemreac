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
- Python (2.7 or >=3.3)
    
In addition to python, the following python packages are required
(versions indicate what is tested):

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

Building and installing
-----------------------
Once non-python prerequisites are in installed, you may procede e.g. as:

::

    $ git clone https://bitbucket.org/bjodah/chemreac.git
    $ cd chemreac
    $ pip install --user --upgrade -r requirements.txt
    $ python setup.py install --user
    $ py.test


the above procedure works on Ubuntu 14.04 for example. See the `Python docs <https://docs.python.org/2/install/index.html#install-index>`_ for more information on how to install e.g. system wide.

To specify an alternative LAPACK lib, set the environment variable LLAPACK, e.g.:

::

    $ LLAPACK=openblas python setup.py build_ext --inplace


Tests
-----
Run ``py.test``, possibly with explicit ``PYTHONPATH`` (e.g. if ``build_ext --inplace`` was used)

::

    $ PYTHONPATH=`pwd`:$PYTHONPATH py.test

All tests should pass (or xfail). If they are not, please `file an issue <https://github.com/bjodah/chemreac/issues>`_.

.. install-end

Status
======


Continuous integration
----------------------
.. ci-start

In order to minimize the risk of (re)introducing bugs the code base, 
it is continuously built on two CI services:

- `travis-ci.org <https://travis-ci.org/bjodah/chemreac>`_
- `drone.io <https://drone.io/github.com/bjodah/chemreac>`_

.. image:: https://travis-ci.org/bjodah/chemreac.png?branch=master
   :target: https://travis-ci.org/bjodah/chemreac

above you find build status shield for travis-ci (Py 2.7, Py 3.4, no OpenMP, runs coveralls).


.. image:: https://drone.io/github.com/bjodah/chemreac/status.png
   :target: https://drone.io/github.com/bjodah/chemreac/latest

above you find build status shield for drone.io (Py 2.7, uses OpenMP, tests sundials backend, 
runs slow tests, build docs as artifact, html coverage report as artifact):


.. ci-end

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
