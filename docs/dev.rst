Development
===========

As an open-source project other researchers/students are encouraged to improve upon
the code/documentation. The code is released under the permissive `BSD 2-Clause
License <http://opensource.org/licenses/BSD-2-Clause>`_. The source
code repository is found at https://github.com/bjodah/chemreac which is also
where issues and patches ("pull requests") are posted.

If you are new to github and the like, a good start is SymPy's wiki
entry on `development workflow
<https://github.com/sympy/sympy/wiki/Development-workflow>`_ (minus
the SymPy specific parts obviously).

Environment
-----------
The documentation assumes a \*NIX like environment, or at
least familiarity with this environment for tweaking build scripts
etc. Feedback from users and developers on other platforms are most
welcome.

Unit tests
----------
The correctness of the implementation is verified through the use of
unit tests. Ideally all code should be covered with tests.

Continuous integration
----------------------
.. include:: ../README.rst
    :start-after: .. ci-start
    :end-before: .. ci-end

Test coverage
-------------

.. image:: https://coveralls.io/repos/bjodah/chemreac/badge.png?branch=master
   :target: https://coveralls.io/r/bjodah/chemreac?branch=master
   :alt: Test coverage

unit tests are a great tool provided they actually cover the majority
of the code base, in order to keep track of the coverage coveralls is
used:


Coding standards
----------------
The Python code should comply with `PEP8
<http://legacy.python.org/dev/peps/pep-0008/>`_. It is also checked by
the CI.
For the C++ code, `PyNE's C/C++ style guide <http://pyne.io/devsguide/style_guide.html#c-c-style-guide>`_ is recommended.

Documentation standard
----------------------
The API docs are generated automatically and `numpydoc
<https://pypi.python.org/pypi/numpydoc>`_ is used to 
automatically parse docstrings of functions and classes. 
See the `NumPy/SciPy documentation guide
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ for more info.
