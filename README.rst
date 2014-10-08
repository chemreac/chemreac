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

the documentation is hostad at

- http://chemreac.readthedocs.org/
- http://bjodah.github.io/chemreac/docs/

Installation
============
.. install-start

Easiest way to install chemreac (on linux) is by using conda to pull it from
https://binstar.org/bjodah/chemreac

Below you will find instructions for installation by building it from source.
You may also look in ``scripts/`` folder for automated install scripts used
by the continuous integration servers.

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
- numpy>=1.9.0
- cython>=0.21.0
- mako>=1.0.0
- quantities>=0.10.1
- pytest>=2.5.2
- pytest-pep8>=1.0.6
- scipy>=0.14.0
- matplotlib>=1.4.0
- periodictable>=1.4.1
- future>=0.13.0
- pycompilation>=0.3.6
- pycodeexport>=0.0.5
- https://github.com/sympy/sympy/archive/master.zip (will req. 0.7.6)

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

In order to minimize the risk of (re)introducing bugs into the code
base, it is continuously built on two CI services:

- `travis-ci.org <https://travis-ci.org/bjodah/chemreac>`_
- `drone.io <https://drone.io/github.com/bjodah/chemreac>`_

.. image:: https://travis-ci.org/bjodah/chemreac.png?branch=master
   :target: https://travis-ci.org/bjodah/chemreac

above you can find the build status shield for travis-ci (Py 2.7, Py
3.4, no OpenMP, runs coveralls, builds docs and pushes them to the
gh-pages branch).


.. image:: https://drone.io/github.com/bjodah/chemreac/status.png
   :target: https://drone.io/github.com/bjodah/chemreac/latest

above you can find build status shield for drone.io (Py 2.7, uses OpenMP, tests sundials backend, build docs as artifact, html coverage report as artifact):


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
