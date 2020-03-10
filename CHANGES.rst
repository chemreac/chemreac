v0.9.1
======
- Update setup.py

v0.9.0
======
- Require sundials>=5.1.0, pycvodes>=0.12.0

v0.8.4
======
- Update pycvodes

v0.8.3
======
- Warn if fields given & doserate/density in vars

v0.8.2
======
- Fixes to ReactionDiffusion.from_ReactionSystem

v0.8.1
======
- Support for constraints (sundials>=0.3.2)

v0.8.0
======
- New AnyODE version (19)
- New upstream pycvodes version (0.11.1)

v0.7.13
=======
- Added support for constraints (sundials>=3.2.0)

v0.7.12
=======
- Updated requirements

v0.7.11
=======
- Back-ported features from v0.8.0

v0.7.10
=======
- Specify (upper) version requirements

v0.7.9
======
- Updated to AnyODE 15 & block_diag_ilu-0.3.6
- Only require C++11

v0.7.8
======
- Added chemreac.util.anim

v0.7.7
======
- AnyODE v14

v0.7.6
======
- More robust build (conda-recipe)
- Update AnyODE

v0.7.5
======
- AnyODE v12
- chemreac._odesys.ODESys now accepts atol in the form of a dict

v0.7.4
======
- New attr for ReactionDiffusion: param_names
- New (provisional) module: ``chemreac._odesys``.

v0.7.3
======
- Use pycvodes-0.9.0 (ABI change)

v0.7.2
======
- Setup fix: explicit package_data in setup.py (conda related)

v0.7.1
======
- Bumpy AnyODE version

v0.7.0
======
- Refactor to use new versions of chempy, AnyODE
- Refactored to use separate block_diag_ilu, finitediff, pycvodes
- Support for base 2 in log-transform.

v0.6.0
======
- ``chemreac.integrate.run`` used solely for kwargs from environment variable:
      ``CHEMREAC_INTEGRATION_KWARGS`` (was ``CHEMREAC_SOLVER_KWARGS``)
- Units handled differently in ``ReactionDiffusion``:
    - ``ReactionDiffusion.__init__`` accepts unitless numbers and a ``unit_registry``
    - ``ReactionDiffusion.nondimensionalisation()`` accepts numbers wiht units and a ``unit_registry``
    - ``ReactionDiffusion.with_units`` accessor.
- ``Integration.nondimensionalistion()`` analogous above.
- Uses ChemPy > 0.4.1
- ``ReactionDiffusion`` can now be used with ``pickle``
- New serialization format (json and pickle)

v0.5.0
======
- neval_f, neval_j -> nfev, njev
- new upstream: block_diag_ilu, pycvodes
- New parameters: ilu_limit, n_jac_diags, first_step
- Change solver choice "sundials" to "cvode" (prepare for arkode support).
- For sundials backend: Don't use iterative 0, 1, 2, 3. Instead use:
   - linear_solver: {'default', 'dense', 'banded', 'gmres', 'bicgstab', 'tfqmr'}
   - iter_type: {'functional', 'newton'}
- Refactored C++ code generation for jacobian routines

v0.4.0
======
- Don't use consants FLAT, CYLINDRICAL, SPHERICAL. Instead use 'f', 'c', 's'
- Drop constants GEOM_ORDER, DENSE, BANDED, SPARSE, GEOMS


v0.3
====
- ReactionSystem.from_ReactionDiffusion, ReactionSystem.to_ReactionDiffusion ->
      ReactionDiffusion.from_ReactionSystem, ReactionDiffusion.to_ReactionSystem
- Use chempy for Substance, Reaction, ReactionSystem etc.
   - ReactionSystem got a new signature
- bin_k_factor/bin_k_factor_span replaced with fields/g_values and modulated_rxns/modulation
- Added support for units (.units, .constants), new serialization format.
- Moved repository to github.com/chemreac/chemreac
- Enabled use of Sundial's iterative linear solvers.
- Added Incomplete LU preconditioner
- Bug fixes (compressed_jac_cmaj)
- Introspection of jacobian/preconditioner (coloured_spy, data dumping)
- New Docker images for CI
- Removed unused code found by pyflakes
- Updated logo and webpage
- Updated setup.py for better automated releases.
- Drop support for k_err and D_err attributes in ReactionDiffusion
- Add support for longtable in chemreac.util.table.rsys2pdf_table(...)
- Clarify if input data to plotting routines should be transformed or not.
- Enable support for file objects in load/dump in serialization
- Improvements to rsys2pdf_table (delete -> save, units)
- Better line color / style selection in plotting functions (28 unique)
- "order" attribute added to chemreac.chemistry.Reaction
