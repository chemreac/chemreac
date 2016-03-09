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
