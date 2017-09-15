Development
===========

As an open-source project, other researchers/students are encouraged to improve upon
the code/documentation. The code is released under the permissive `BSD 2-Clause
License <http://opensource.org/licenses/BSD-2-Clause>`_. The source
code repository is found at https://github.com/chemreac/chemreac which is also
where issues and patches ("pull requests") are accepted.


Coding standards
----------------
The Python code should comply with `PEP8
<http://legacy.python.org/dev/peps/pep-0008/>`_.  Before submitting a
Pull Request you also want to pass ``--pep8`` to ``py.test`` (it is
done automatically by the ``./scripts/run_tests.sh`` script but you
need ``pytest-pep8`` for it to work).

PEP8 complience is
also checked by the CI servers. For the C++ code, `PyNE's C/C++ style
guide <http://pyne.io/devsguide/style_guide.html#c-c-style-guide>`_ is
recommended.

Documentation standard
----------------------
The API docs are generated automatically and `numpydoc
<https://pypi.python.org/pypi/numpydoc>`_ is used to
automatically parse docstrings of functions and classes.
See the `NumPy/SciPy documentation guide
<https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_ for more info.

Unit tests
----------
The correctness of the implementation is verified through the use of
(hundreds of) unit tests. Ideally all code should be covered with tests.

Continuous integration
~~~~~~~~~~~~~~~~~~~~~~
.. include:: ../README.rst
    :start-after: .. ci-start
    :end-before: .. ci-end

Test coverage
~~~~~~~~~~~~~
Unit testing is a great tool, provided the tests actually cover the majority
of the code base. In order to keep track of the test coverage a
`report <http://hera.physchem.kth.se/~chemreac/branches/master/htmlcov>`_ is
generated for each change of the code.

Environment
-----------
In theory ``chemreac`` should be cross-platform. However, all development (and
testing) is done on Linux, and the documentation assumes a posix compliant
environment, or at least familiarity with this environment. Feedback from users
and developers on other platforms is most welcome.

If you are new to `github <https://github.com/>`_ and the like, a good start
is SymPy's wiki entry on `development workflow
<https://github.com/sympy/sympy/wiki/Development-workflow>`_ (minus the SymPy
specific parts obviously).

Example of setting up a development environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here is an example to get you started developing chemreac and contributing
to the code:

1. `Fork <https://help.github.com/articles/fork-a-repo/>`_ the repository on github
2. Install the dependencies specified in the README.rst file
3. Clone your own repo to get a local copy:

::

    $ git clone https://github.com/YOUR_GITHUB_USERNAME/chemreac.git

4. Optionally set-up the pre-commit hook (prevents commiting untested code).
::

    $ cd .git/hooks
    $ ln -s ../../scripts/pre-commit.sh pre-commit
    $ cd -

5. Create a new branch
::

    $ git checkout -b fix_silly_bug


6. Edit files, add test(s), make sure all tests pass by running:
::

    $ ./scripts/run_tests.sh


7. make a commit, and push your changes:
::

    $ git commit -am "Fixed a silly bug in the complex-thingy."
    $ git push origin --set-upstream fix_silly_bug


8. Go to your forked repo on github and create a `pull-request <https://help.github.com/articles/using-pull-requests/>`_ from there.
