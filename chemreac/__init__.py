# -*- coding: utf-8 -*-
"""
Python extension for reaction diffusion.
"""

import os
import subprocess

#     Debugging:
#
# ===== testing package: chemreac-0.2.2.dev-np19py27_0 =====
# Importing chemreac...
# linux-vdso.so.1 =>  (0x00007fffbedfe000)
# libsundials_cvode.so.1 => /usr/local/lib/libsundials_cvode.so.1 (0x00
# 007f95362b8000)
# libopenblas.so.0 => /usr/local/lib/libopenblas.so.0 (0x00007f95353f00
# 00)
# libsundials_nvecserial.so.0 => /usr/local/lib/libsundials_nvecserial.
# so.0 (0x00007f95351eb000)
# libpython2.7.so.1.0 => /home/bjorn/miniconda/envs/_test/lib/python2.7
# /site-packages/chemreac/../../../libpython2.7.so.1.0 (0x00007f9534e09000)
# libgfortran.so.3 => /usr/lib/x86_64-linux-gnu/libgfortran.so.3 (0x000
# 07f9534ac8000)
# libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007f9
# 5347c4000)
# libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007f95345ae
# 000)
# libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f95341e7000)
# libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f9533ee1000)
# libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007f953
# 3cc3000)
# libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007f9533abe000)
# libutil.so.1 => /lib/x86_64-linux-gnu/libutil.so.1 (0x00007f95338bb00
# 0)
# libquadmath.so.0 => /usr/lib/x86_64-linux-gnu/libquadmath.so.0 (0x000
# 07f953367e000)
# /lib64/ld-linux-x86-64.so.2 (0x00007f9536732000)

# Importing chemreac failed
# /home/bjorn/miniconda/envs/_test/lib/python2.7/site-packages/chemreac
# /_chemreac.so: undefined symbol: __exp_finite
# Traceback (most recent call last):
#   File "/home/bjorn/miniconda/conda-bld/test-tmp_dir/run_test.py", li
# ne 34, in <module>
#     import chemreac
#   File "/home/bjorn/miniconda/envs/_test/lib/python2.7/site-packages/
# chemreac/__init__.py", line 15, in <module>
#     from .core import (
#   File "/home/bjorn/miniconda/envs/_test/lib/python2.7/site-packages/
# chemreac/core.py", line 17, in <module>
#     from ._chemreac import CppReactionDiffusion
# ImportError: /home/bjorn/miniconda/envs/_test/lib/python2.7/site-pack
# ages/chemreac/_chemreac.so: undefined symbol: __exp_finite
# Traceback (most recent call last):
#   File "/home/bjorn/miniconda/conda-bld/test-tmp_dir/run_test.py", li
# ne 39, in <module>
#     raise e
# ImportError: /home/bjorn/miniconda/envs/_test/lib/python2.7/site-pack
# ages/chemreac/_chemreac.so: undefined symbol: __exp_finite
# TESTS FAILED: chemreac-0.2.2.dev-np19py27_0

import glob
_path = os.path.join(
    os.path.dirname(__file__),
    '_chemreac*.so'
)
print(subprocess.check_output(['ldd', '-r', glob.glob(_path)[0]]))

from .release import __version__

from .core import (
    ReactionDiffusion, DENSE, BANDED, SPARSE, FLAT,
    CYLINDRICAL, SPHERICAL, Geom_names
)
