========
chemreac
========

.. image:: https://drone.io/bitbucket.org/bjodah/chemreac/status.png
   :target: https://drone.io/bitbucket.org/bjodah/chemreac/latest

chremreac aims to be a library collecting tools and utilities for modeling 
of chemical kinetics problems, primarily in the context of aqueous phase 
with external radiation fields.

the documentation is found at [INSERT URL HERE...]

.. install-start

Installation
============
chemreac requires:

    C++ compiler with full C++11 support (e.g. GCC >= 4.8)
    Fortran compiler with ISO_C_BINDING support (Fortran 2003 standard) (e.g. gfortran)
    LAPACK
    Sundials 2.5
    Python_ (2.7)
    
In addition to python, the following python packages are requiered:

   #. argh>=0.25.0
   #. numpy>=1.8.1
   #. cython>=0.19.1
   #. mako>=0.5.0
   #. quantities>=0.10.1
   #. pytest>=2.5.2
   #. scipy>=0.9.0
   #. periodictable>=1.4.1
   #. future>=0.12.3
   #. https://github.com/bjodah/pycompilation/archive/v0.2.21.tar.gz
   #. https://github.com/sympy/sympy/archive/master.zip
   #. `NumPy <http://www.numpy.org/>`_ (>= 1.8.0)
   #. `SciPy <http://www.scipy.org/>`_
   #. `Cython <http://cython.org/>`_ (>= 0.19.1)

is needed. For rendering the documentation you also need:

   #. `Sphinx <http://sphinx-doc.org/>`_
      

to run all the tests you also need:

    valgrind
    graphviz
    dot2tex

    git clone https://bitbucket.org/bjodah/chemreac.git
    cd chemreac
    pip install --user --upgrade -r requirements.txt
    python setup.py build_ext --inplace
    py.test

should be enough (tested on Ubuntu 12.04). Add to $PYTHONPATH to use
or ``python setup.py install --user``

To specify an alternative LAPACK lib set the environment variable LLAPACK, e.g.:

```
    LLAPACK=openblas python setup.py build_ext --inplace
```

You may also look at ``scrips/test_drone.io.sh`` for the CI install.

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

- C++11 compliant compiler

.. install-end

TODO
======
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
Björn Dahlgren, contact (gmail adress): bjodah
