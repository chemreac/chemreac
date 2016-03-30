========
chemreac
========

.. image:: http://hera.physchem.kth.se:9090/api/badges/chemreac/chemreac/status.svg
   :target: http://hera.physchem.kth.se:9090/chemreac/chemreac
   :alt: Build status
.. image:: https://img.shields.io/pypi/v/chemreac.svg
   :target: https://pypi.python.org/pypi/chemreac
   :alt: PyPI version
.. image:: https://img.shields.io/badge/python-2.7,3.4,3.5-blue.svg
   :target: https://www.python.org/
   :alt: Python version
.. image:: https://zenodo.org/badge/8840/chemreac/chemreac.svg
   :target: https://zenodo.org/badge/latestdoi/8840/chemreac/chemreac
   :alt: DOI (digital object identifier) for latest release on Zenodo
.. image:: https://img.shields.io/pypi/l/chemreac.svg
   :target: https://github.com/bjodah/chemreac/blob/master/LICENSE
   :alt: License
.. image:: http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat
   :target: http://hera.physchem.kth.se/~chemreac/benchmarks
   :alt: airspeedvelocity
.. image:: http://hera.physchem.kth.se/~chemreac/branches/master/htmlcov/coverage.svg
   :target: http://hera.physchem.kth.se/~chemreac/branches/master/htmlcov
   :alt: coverage

chremreac is an open source library which aims to collect tools and utilities for
modeling of chemical kinetics problems. It is primarily designed to
be useful in the context of aqueous phase with external radiation fields.

The (heavy) numerical processing is done in rountines written in C++ which have
been wrapped in a `Python <https://www.python.org>`_ interface using
`Cython <https://www.cython.org>`_.

Documentation
=============

- Documentation for the latest stable release is found here:
`<https://pythonhosted.org/chemreac>`_
- and development documentation for the current master branch is found here:
`<http://hera.physchem.kth.se/~chemreac/branches/master/html>`_).


Installation
============
.. install-start

The easiest way to install ``chemreac`` (on linux) is to use
`conda <http://docs.continuum.io/anaconda/index.html>`_:

::

   $ conda install -c chemreac chemreac pytest

and you're done! To check if ``chemreac`` is installed correctly you may run:

::

    $ pytest --pyargs chemreac

which should run the test suite (all tests should pass or xfail).
If you are using a special platform with non-standard math libraries you
may need to compile your own pacakges. The "recipes" for the conda packages
are kept `here <https://github.com/chemreac/chemreac_anaconda/>`_.

If you are not using the conda package manager you can still install
``chemreac`` from source. You will find the instructions for doing so below.

Building from source
--------------------
Below you will find instructions for installation by building ``chemreac`` from source.
You may also look in ``scripts/`` folder for automated install scripts used
on the continuous integration servers.

Prerequisites
~~~~~~~~~~~~~
Version numbers of dependencies indicate what has been tested:

- C++ compiler with C++11 support (e.g. `GCC <https://gcc.gnu.org/>`_ >= 4.8)
- LAPACK (provided by e.g. `OpenBLAS <http://www.openblas.net/>`_)
- `Sundials <http://computation.llnl.gov/casc/sundials/main.html>`_ 2.6.2
- `Python <https://www.python.org>`_ (2.7 or >=3.4)

In addition to the standard library provided by Python, a number of python
packages are required (see `setup.py <./setup.py>`_), they are found on `PyPI
<https://pypi.python.org/pypi>`_ and are automatically installed when
using ``pip``.

For rendering the documentation you also need:

- `Sphinx <http://sphinx-doc.org/>`_
- `numpydoc <https://pypi.python.org/pypi/numpydoc>`_
- `sphinx_rtd_theme <https://pypi.python.org/pypi/sphinx_rtd_theme>`_

to run all the tests you also need these tools:

- valgrind
- graphviz
- dot2tex
- pdflatex


Building and installing
~~~~~~~~~~~~~~~~~~~~~~~
Once non-python prerequisites are in installed, you may procede e.g. as:

::

    $ git clone https://github.com/chemreac/chemreac.git
    $ cd chemreac
    $ pip install --user -e .[all]
    $ ./scripts/run_tests.sh


the above procedure works on Ubuntu 14.04 for example.

To specify an alternative LAPACK lib, set the environment variable LLAPACK, e.g.:

::

    $ LLAPACK=openblas python setup.py build_ext --inplace

The following environment variables are also supported by
``setup.py`` (defaults to "0", enable by setting them to "1"):

+-----------------------------------------+-------+--------------------------------------------------+
|Environament variable                    |Default|Action                                            |
+=========================================+=======+==================================================+
|``WITH_OPENMP``                          |0      |Enable parallell evaluation of rhs and jac.       |
+-----------------------------------------+-------+--------------------------------------------------+
|``WITH_BLOCK_DIAG_ILU_DGETRF``           |0      |Use unblocked version of dgetrf instead of LAPACK |
+-----------------------------------------+-------+--------------------------------------------------+
|``WITH_BLOCK_DIAG_ILU_OPENMP``           |0      |Evaluate LU decomposition of blocks in parallel   |
+-----------------------------------------+-------+--------------------------------------------------+
|``WITH_DATA_DUMPING``                    |0      |For debugging purposes only                       |
+-----------------------------------------+-------+--------------------------------------------------+
|``WITH_DEBUG``                           |0      |For debugging purposes only                       |
+-----------------------------------------+-------+--------------------------------------------------+

Enabling the first three is known to provide significant speed up for some scenarios (performance is
system dependent, hence recommendations are not possible to give without benchmarking).

Tests
-----
If you have ``py.test`` installed you may run the test suite on the
installed package:

::

    $ py.test --pyargs chemreac

All tests should pass (or xfail). If they do not, please `file an
issue <https://github.com/chemreac/chemreac/issues>`_.

.. install-end

Status
======
Both the correctness (continuous integration) and the performance
(benchmarks) of the code base is monitored continuously.

Continuous integration
----------------------
.. ci-start

In order to minimize the risk of (re)introducing bugs into the code
base, it is continuously built on a CI server:

.. image:: http://hera.physchem.kth.se:9090/api/badges/chemreac/chemreac/status.svg
   :target: http://hera.physchem.kth.se:9090/chemreac/chemreac
   :alt: Build status

.. ci-end

Performace tracking
-------------------
Benchmarks for tracking the performance of the library are kept at
https://github.com/chemreac/chemreac-benchmarks


License
=======
The source code is Open Source and is released under the very permissive
"simplified (2-clause) BSD license". See ``LICENSE.txt`` for further details.
Contributors are welcome to suggest improvements at https://github.com/chemreac/chemreac

Author
======
Björn Dahlgren, contact:
 - gmail address: bjodah
 - kth.se address: bda
