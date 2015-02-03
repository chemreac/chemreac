#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import sys

from distutils.core import setup, Command

pkg_name = 'chemreac'

with open(os.path.join(pkg_name,'__init__.py')) as f:
    long_description = f.read().split('"""')[1]

# Reading the version is a bit tricky: the same commit could actually
# correspond to multiple versions. e.g. a commit could be tagged:
# v0.1.0-rc1, v0.1.0-rc2, v0.1.0
# Hence, the build environment needs a way to override the version
# string found by setup.py (unless setup.py is to depend on git)
# the solution right now is for setup.py
# to look for the environment variable CHEMREAC_RELEASE_VERSION
# if not set it will exec the contents of: ./chemreac/release.py
CHEMREAC_RELEASE_VERSION = os.environ.get('CHEMREAC_RELEASE_VERSION', '')

release_py_path = os.path.join(pkg_name, 'release.py')

if len(CHEMREAC_RELEASE_VERSION) > 1:
    __version__ = CHEMREAC_RELEASE_VERSION
else:
    # read __version__ attribute:
    exec(open(release_py_path).read())

with open(pkg_name+'/__init__.py') as f:
    long_description = f.read().split('"""')[1]

DEBUG = True if os.environ.get('USE_DEBUG', False) else False
USE_OPENMP = True if os.environ.get('USE_OPENMP', False) else False
LLAPACK = os.environ.get('LLAPACK', 'lapack')
WITH_BLOCK_DIAG_ILU_DGETRF = os.environ.get('WITH_BLOCK_DIAG_ILU_DGETRF', '0') == '1'
WITH_BLOCK_DIAG_ILU_OPENMP = os.environ.get('WITH_BLOCK_DIAG_ILU_OPENMP', '0') == '1'
WITH_DATA_DUMPING = os.environ.get('WITH_DATA_DUMPING', '0') == '1'

CONDA_BUILD = os.environ.get('CONDA_BUILD', '0') == '1'
ON_RTD = os.environ.get('READTHEDOCS', None) == 'True'
ON_DRONE = os.environ.get('DRONE', 'false') == 'true'
ON_TRAVIS = os.environ.get('TRAVIS', 'flse') == 'true'

if CONDA_BUILD:
    open('__conda_version__.txt', 'w').write(__version__)

flags = []
options = ['pic', 'warn']
if DEBUG:
    print("Building chemreac with debugging enabled.")
    options += ['debug']
else:
    if not (ON_DRONE or ON_TRAVIS):
        if CONDA_BUILD:
            # -ffast-math buggy in anaconda
            flags += ['-O2', '-funroll-loops']
        else:
            options += ['fast']  # -ffast-math -funroll-loops

if WITH_BLOCK_DIAG_ILU_OPENMP:
    options += ['openmp']

cmdclass_ = {}

IDEMPOTENT_INVOCATION = False
if len(sys.argv) > 1:
    if '--help' in sys.argv[1:] or sys.argv[1] in (
            '--help-commands', 'egg_info', 'clean', '--version'):
        IDEMPOTENT_INVOCATION = True
elif len(sys.argv) == 1:
    IDEMPOTENT_INVOCATION = True

if ON_RTD or IDEMPOTENT_INVOCATION:
    # Enbale pip to probe setup.py before all requirements are installed
    ext_modules_ = []
else:
    import pickle
    import numpy as np
    template_path = 'src/chemreac_template.cpp'
    rendered_path = 'src/chemreac.cpp'
    USE_TEMPLATE = not os.path.exists(rendered_path)
    try:
        from pycodeexport.dist import PCEExtension, pce_build_ext, pce_sdist
    except ImportError:
        if USE_TEMPLATE:
            print("This is not source distribution. pycodeexport is needed:", sys.exc_info()[0])
            raise
        # If building from sdist no need for more than pycompilation
        from pycompilation.dist import PCExtension as PCEExtension
        from pycompilation.dist import pc_build_ext as pce_build_ext
        from pycompilation.dist import pc_sdist as pce_sdist

    cmdclass_['build_ext'] = pce_build_ext
    cmdclass_['sdist'] = pce_sdist
    subsd = {'USE_OPENMP': USE_OPENMP}
    pyx_path = 'chemreac/_chemreac.pyx'
    using_pyx = os.path.exists(pyx_path)
    pyx_or_cpp = pyx_path if using_pyx else pyx_path[:-3]+'cpp'
    sources = [
        template_path if USE_TEMPLATE else rendered_path,
        'src/finitediff/finitediff/fornberg.f90',
        'src/finitediff/finitediff/c_fornberg.f90',
        pyx_or_cpp,
    ]

    ext_modules_ = [
        PCEExtension(
            "chemreac._chemreac",
            sources=sources,
            template_regexps=[
                (r'^(\w+)_template.(\w+)$', r'\1.\2', subsd),
            ] if USE_TEMPLATE else (),
            pycompilation_compile_kwargs={
                'per_file_kwargs': {
                    'src/chemreac.cpp': {
                        'std': 'c++0x',
                        # 'fast' doesn't work on drone.io
                        'flags': flags,
                        'options': options +
                        (['openmp'] if USE_OPENMP else []),
                        'define': [] +
                        (['DEBUG'] if DEBUG else []) +
                        (['WITH_DATA_DUMPING'] if
                         WITH_DATA_DUMPING else []) +
                        (['WITH_BLOCK_DIAG_ILU_OPENMP'] if
                         WITH_BLOCK_DIAG_ILU_OPENMP else []) +
                        (['WITH_BLOCK_DIAG_ILU_DGETRF'] if
                         WITH_BLOCK_DIAG_ILU_DGETRF else []),
                    },
                    'src/chemreac_sundials.cpp': {
                        'std': 'c++0x',
                        'flags': flags,
                        'options': options
                    },
                    pyx_or_cpp: {
                        'cy_kwargs': {'annotate': True},
                        'std': 'c++0x',
                        'gdb_debug': DEBUG
                    } if using_pyx else {
                        'std': 'c++0x',
                    }
                },
                'flags': flags,
                'options': options,
            },
            pycompilation_link_kwargs={
                'options': (['openmp'] if USE_OPENMP else []),
                'std': 'c++0x',
            },
            include_dirs=['src/', 'src/finitediff/finitediff/',
                          np.get_include()],
            libraries=['sundials_cvodes', LLAPACK, 'sundials_nvecserial', 'm'],
            logger=True,
        )
    ]

modules = [
    pkg_name+'.util',
]

tests = [
    pkg_name+'.tests',
    pkg_name+'.util.tests',
]

package_data = {
    pkg_name: ['tests/*.json', 'tests/*.txt']
}

classifiers = [
    "Development Status :: 3 - Alpha",
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
]

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description='Python package for modeling chemical kinetics with diffusion and drift.',
    long_description=long_description,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://github.com/chemreac/' + pkg_name,
    packages=[pkg_name] + modules + tests,
    package_data=package_data,
    cmdclass=cmdclass_,
    ext_modules=ext_modules_,
)

if __name__ == '__main__':
    try:
        if len(CHEMREAC_RELEASE_VERSION) > 1:
            # Same commit should generate different sdist
            # depending on tagged version (set CHEMREAC_RELEASE_VERSION)
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write("__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if len(CHEMREAC_RELEASE_VERSION) > 1:
            shutil.move(release_py_path+'__temp__', release_py_path)
