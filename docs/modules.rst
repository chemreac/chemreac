Modules Reference
=================

At package level names that are gauranteed to be exported are: ``__version__`` and ``ReactionDiffusion``. For example use of the different modules you may look at the examples.

There are a number of environment variables which may be set to change the defaults
of some runtime parameters. Here is an overview:

CHEMREAC_INTEGRATOR: sundials, e.g. 'scipy'
CHEMREAC_INTEGRATOR_KWARGS: nil, e.g. '{"iterative": 1}'


.. toctree::
   :maxdepth: 4

   core.rst
   integrate.rst
   chemistry.rst
   util/index.rst
