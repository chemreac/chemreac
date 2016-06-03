#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import shutil
import sys

from setuptools import setup

pkg_name = 'chemreac'


def _path_under_setup(*args):
    return os.path.join(os.path.dirname(__file__), *args)

# Reading the version is a bit tricky: the same commit could actually
# correspond to multiple versions. e.g. a commit could be tagged both:
# v0.1.0-rc1, v0.1.0
#
# Hence, the build environment needs a way to determine the version based
# on meta data not contained in the version controlled files.
# This could be done either by having setup.py inspect git tags
# (which makes setup.py dependent on git), or have setup.py look for an
# environment variable.
#
# The latter method is used for now:
# The variable $CHEMREAC_RELEASE_VERSION is matched against "v*" and
# if valid __version__ is set accordingly. If mathcing fails setup.py
# will exec the contents of: ./chemreac/_release.py
#
# To complicate things further conda-build drops most environment
# variables, so for conda based builds to work need setup.py to write
# the version as a string to a file named '__conda_version__.txt'

RELEASE_VERSION = os.environ.get('CHEMREAC_RELEASE_VERSION', '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
CONDA_BUILD = os.environ.get('CONDA_BUILD', '0') == '1'
if CONDA_BUILD:
    try:
        RELEASE_VERSION = 'v' + io.open('__conda_version__.txt', 'rt',
                                        encoding='utf-8').readline().rstrip()
    except IOError:
        pass

release_py_path = _path_under_setup(pkg_name, '_release.py')

if len(RELEASE_VERSION) > 1 and RELEASE_VERSION[0] == 'v':
    TAGGED_RELEASE = True
    __version__ = RELEASE_VERSION[1:]
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(io.open(release_py_path, encoding='utf-8').read())

WITH_OPENMP = os.environ.get('WITH_OPENMP', '0') == '1'
LLAPACK = os.environ.get('LLAPACK', 'lapack')
WITH_BLOCK_DIAG_ILU_DGETRF = os.environ.get(
    'WITH_BLOCK_DIAG_ILU_DGETRF', '0') == '1'
WITH_BLOCK_DIAG_ILU_OPENMP = os.environ.get(
    'WITH_BLOCK_DIAG_ILU_OPENMP', '0') == '1'
WITH_DATA_DUMPING = os.environ.get('WITH_DATA_DUMPING', '0') == '1'
WITH_DEBUG = os.environ.get('WITH_DEBUG', '0') == '1'

ON_DRONE = os.environ.get('DRONE', 'false') == 'true'
ON_TRAVIS = os.environ.get('TRAVIS', 'flse') == 'true'

# See pycompilation for details on "options"
options = ['pic', 'warn']
if WITH_DEBUG:
    print("Building chemreac with debugging enabled.")
    options += ['debug']
    flags = []
else:
    flags = ['-O3']
    if not (ON_DRONE or ON_TRAVIS):
        if CONDA_BUILD:
            # -ffast-math buggy in anaconda
            flags += ['-funroll-loops']
        # else:
        #     options += ['fast']  # -ffast-math -funroll-loops

if os.environ.get('WITH_PROFILE', '0') == '1':
    flags += ['-g', '-pg']

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

if IDEMPOTENT_INVOCATION:
    # Enbale pip to probe setup.py before all requirements are installed
    ext_modules_ = []
else:
    import pickle
    import numpy as np
    template_path = 'src/chemreac_template.cpp'
    rendered_path = 'src/chemreac.cpp'
    # Source distributions contain rendered sources
    USE_TEMPLATE = os.path.exists(template_path)
    try:
        from pycodeexport.dist import PCEExtension, pce_build_ext, pce_sdist
    except ImportError:
        if USE_TEMPLATE:
            print("This is not source distribution. pycodeexport is needed:",
                  sys.exc_info()[0])
            raise
        # If building from sdist no need for more than pycompilation
        from pycompilation.dist import PCExtension as PCEExtension
        from pycompilation.dist import pc_build_ext as pce_build_ext
        from pycompilation.dist import pc_sdist as pce_sdist

    cmdclass_['build_ext'] = pce_build_ext
    cmdclass_['sdist'] = pce_sdist
    subsd = {'WITH_OPENMP': WITH_OPENMP}
    pyx_path = 'chemreac/_chemreac.pyx'
    using_pyx = os.path.exists(pyx_path)  # No pyx in source dist
    pyx_or_cpp = pyx_path if using_pyx else pyx_path[:-3]+'cpp'
    sources = [
        template_path if USE_TEMPLATE else rendered_path,
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
                        'std': 'c++11',
                        # 'fast' doesn't work on drone.io
                        'flags': flags,
                        'options': options +
                        (['openmp'] if WITH_OPENMP else []),
                        'define': [] +
                        (['WITH_DEBUG'] if WITH_DEBUG else []) +
                        (['WITH_DATA_DUMPING'] if
                         WITH_DATA_DUMPING else []) +
                        (['WITH_BLOCK_DIAG_ILU_OPENMP'] if
                         WITH_BLOCK_DIAG_ILU_OPENMP else []) +
                        (['WITH_BLOCK_DIAG_ILU_DGETRF'] if
                         WITH_BLOCK_DIAG_ILU_DGETRF else []),
                    },
                    'src/chemreac_sundials.cpp': {
                        'std': 'c++11',
                        'flags': flags,
                        'options': options
                    },
                    pyx_or_cpp: {
                        'cy_kwargs': {
                            'annotate': True,
                            'embedsignature': True,
                            'linetrace': not TAGGED_RELEASE},
                        'std': 'c++11',
                        'define': ['CYTHON_TRACE=1'],
                        'gdb_debug': WITH_DEBUG
                    } if using_pyx else {
                        'std': 'c++11',
                        'define': ['CYTHON_TRACE=1'],
                        'inc_py': True,
                    }
                },
                'flags': flags,
                'options': options,
            },
            pycompilation_link_kwargs={
                'options': ((['WITH_DEBUG'] if WITH_DEBUG else []) +
                            (['openmp'] if WITH_OPENMP else [])),
                'std': 'c++11',
            },
            include_dirs=['src/', 'src/finitediff/include/',
                          'src/finitediff/external/newton_interval/include/',
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

long_description = io.open('README.rst', encoding='utf-8').read()
with io.open(_path_under_setup(pkg_name, '__init__.py'),
             encoding='utf-8') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
assert len(short_description) > 10 and len(short_description) < 256

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=short_description,
    long_description=long_description,
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    license='BSD',
    keywords=["chemical kinetics", "Smoluchowski equation",
              "advection-diffusion-reaction"],
    url='https://github.com/chemreac/' + pkg_name,
    packages=[pkg_name] + modules + tests,
    package_data=package_data,
    cmdclass=cmdclass_,
    ext_modules=ext_modules_,
    classifiers=classifiers,
    setup_requires=['pycompilation', 'pycodeexport', 'mako'],
    install_requires=['numpy', 'chempy>=0.4.1', 'quantities', 'block_diag_ilu'],
    extras_require={'all': ['argh', 'pytest', 'scipy', 'matplotlib', 'mpld3',
                            'sympy', 'pyodeint', 'pygslodeiv2', 'batemaneq']}

)

if __name__ == '__main__':
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist
            # depending on tagged version (set CHEMREAC_RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write(
                "__version__ = '{}'\n".format(__version__))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__', release_py_path)
