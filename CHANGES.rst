v0.3
====
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
