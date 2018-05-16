Modules Reference
=================

At package level, variables that are gauranteed to be exported are:
``__version__`` and ``ReactionDiffusion``.
For examples of the modules' use, you may look at the examples.

There are a number of environment variables which may be set to change the defaults
of some runtime parameters. Here is an overview:

+--------------------------------+--------------+----------------------+
| Environment variable           | Default      | Example              |
+================================+==============+======================+
| ``CHEMREAC_INTEGRATOR_KWARGS`` | `nil`        | ``{"iterative": 1}`` |
+--------------------------------+--------------+----------------------+

.. toctree::
   :maxdepth: 4

   core.rst
   integrate.rst
   chemistry.rst
   util/index.rst
