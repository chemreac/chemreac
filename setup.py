#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import pprint
import re
import shutil
import subprocess
import sys
import warnings

from setuptools import setup

pkg_name = 'chemreac'
url = 'https://github.com/chemreac/' + pkg_name
license = 'BSD'


def _path_under_setup(*args):
    return os.path.join(os.path.dirname(__file__), *args)


release_py_path = _path_under_setup(pkg_name, '_release.py')
config_py_path = _path_under_setup(pkg_name, '_config.py')
env = None  # silence pyflakes, 'env' is actually set on the next line
exec(open(config_py_path).read())
for k, v in list(env.items()):
    env[k] = os.environ.get('%s_%s' % (pkg_name.upper(), k), v)

_version_env_var = '%s_RELEASE_VERSION' % pkg_name.upper()
RELEASE_VERSION = os.environ.get(_version_env_var, '')

# http://conda.pydata.org/docs/build.html#environment-variables-set-during-the-build-process
CONDA_BUILD = os.environ.get('CONDA_BUILD', '0') == '1'
if CONDA_BUILD:
    try:
        RELEASE_VERSION = 'v' + io.open(
            '__conda_version__.txt', 'rt', encoding='utf-8').readline().rstrip()
    except IOError:
        pass


if len(RELEASE_VERSION) > 1 and RELEASE_VERSION[0] == 'v':
    TAGGED_RELEASE = True
    __version__ = RELEASE_VERSION[1:]
else:
    TAGGED_RELEASE = False
    # read __version__ attribute from _release.py:
    exec(io.open(release_py_path, encoding='utf-8').read())
    if __version__.endswith('git'):
        try:
            _git_version = subprocess.check_output(
                ['git', 'describe', '--dirty']).rstrip().decode('utf-8').replace('-dirty', '.dirty')
        except subprocess.CalledProcessError:
            warnings.warn("A git-archive is being installed - version information incomplete.")
        else:
            if 'develop' not in sys.argv:
                warnings.warn("Using git to derive version: dev-branches may compete.")
                __version__ = re.sub('v([0-9.]+)-(\d+)-(\w+)', r'\1.post\2+\3', _git_version)  # .dev < '' < .post

ON_DRONE = os.environ.get('DRONE', 'false') == 'true'
ON_TRAVIS = os.environ.get('TRAVIS', 'flse') == 'true'

# See pycompilation for details on "options"
options = ['pic', 'warn']
_WITH_DEBUG = env['WITH_DEBUG'] == '1'
_WITH_OPENMP = env['WITH_OPENMP'] == '1'
_WITH_DATA_DUMPING = env['WITH_DATA_DUMPING'] == '1'

if _WITH_DEBUG:
    warnings.warn("Building chemreac with debugging enabled.")
    options += ['debug']
    flags = []
else:
    flags = ['-O3']
    if not (ON_DRONE or ON_TRAVIS):
        if CONDA_BUILD:
            # -ffast-math buggy in anaconda
            flags += ['-funroll-loops']

if _WITH_OPENMP or os.environ.get('BLOCK_DIAG_ILU_WITH_OPENMP', '0') == '1':
    options += ['openmp']

cmdclass_ = {}

# Source distributions contain rendered sources
template_path = 'src/chemreac_template.cpp'
rendered_path = 'src/chemreac.cpp'
USE_TEMPLATE = os.path.exists(template_path)
setup_requires=['pycompilation', 'pycodeexport', 'mako', 'block_diag_ilu', 'pycvodes', 'finitediff'],
install_requires = ['numpy', 'chempy>=0.4.1', 'quantities', 'block_diag_ilu', 'pycvodes', 'finitediff']
package_include = os.path.join(pkg_name, 'include')


if len(sys.argv) > 1 and '--help' not in sys.argv[1:] and sys.argv[1] not in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    import pickle
    import numpy as np
    import finitediff as fd
    import pycvodes as pc
    import block_diag_ilu as bdi
    try:
        from pycodeexport.dist import PCEExtension, pce_build_ext, pce_sdist
    except ImportError:
        if USE_TEMPLATE:
            sys.stderr.write("This is not source distribution. pycodeexport is needed.")
            raise
        # If building from sdist no need for more than pycompilation
        from pycompilation.dist import PCExtension as PCEExtension
        from pycompilation.dist import pc_build_ext as pce_build_ext
        from pycompilation.dist import pc_sdist as pce_sdist

    cmdclass_['build_ext'] = pce_build_ext
    cmdclass_['sdist'] = pce_sdist
    subsd = {'WITH_OPENMP': _WITH_OPENMP}
    pyx_path = 'chemreac/_chemreac.pyx'
    using_pyx = os.path.exists(pyx_path)  # No pyx in source dist
    pyx_or_cpp = pyx_path if using_pyx else pyx_path[:-3]+'cpp'
    sources = [
        template_path if USE_TEMPLATE else rendered_path,
        pyx_or_cpp,
    ]
    _inc_dirs = [
        np.get_include(), fd.get_include(), bdi.get_include(),
        pc.get_include(), package_include, os.path.join('external', 'anyode', 'include')
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
                        'options': options,
                        'define': (['CHEMREAC_WITH_DEBUG'] if _WITH_DEBUG else []) +
                        (['CHEMREAC_WITH_DATA_DUMPING'] if _WITH_DATA_DUMPING else []) +
                        (['BLOCK_DIAG_ILU_WITH_OPENMP'] if os.environ.get('BLOCK_DIAG_ILU_WITH_OPENMP', '') == '1' else []) +
                        (['BLOCK_DIAG_ILU_WITH_DGETRF'] if os.environ.get('BLOCK_DIAG_ILU_WITH_DGETRF', '') == '1' else []),
                    },
                    'src/chemreac_sundials.cpp': {
                        'std': 'c++11',
                        'flags': flags,
                        'options': options
                    },
                    pyx_or_cpp: {
                        'cy_kwargs': {
                            'annotate': True,
                            'include_path': _inc_dirs
                        },
                        'std': 'c++11',
                        'define': ['CYTHON_TRACE=1'],
                        'gdb_debug': _WITH_DEBUG,
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
                'options': options,
                'std': 'c++11',
            },
            include_dirs=_inc_dirs,
            libraries=pc.config['SUNDIALS_LIBS'].split(',') + pc.config['LAPACK'].split(',') + ['m'],
            logger=True,
        )
    ]
else:
    # Enbale pip to probe setup.py before all requirements are installed
    ext_modules_ = []
    if USE_TEMPLATE:
        install_requires += setup_requires

modules = [
    pkg_name+'.util',
]

tests = [
    pkg_name+'.tests',
    pkg_name+'.util.tests',
]

classifiers = [
    "Development Status :: 3 - Alpha",
    'License :: OSI Approved :: BSD License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Mathematics',
]

with io.open(_path_under_setup(pkg_name, '__init__.py'), 'rt', encoding='utf-8') as f:
    short_description = f.read().split('"""')[1].split('\n')[1]
if not 10 < len(short_description) < 255:
    warnings.warn("Short description from __init__.py proably not read correctly.")
long_description = io.open(_path_under_setup('README.rst'),
                           encoding='utf-8').read()
if not len(long_description) > 100:
    warnings.warn("Long description from README.rst probably not read correctly.")
_author, _author_email = io.open(_path_under_setup('AUTHORS'), 'rt', encoding='utf-8').readline().split('<')

setup_kwargs = dict(
    name=pkg_name,
    version=__version__,
    description=short_description,
    long_description=long_description,
    author=_author.strip(),
    author_email=_author_email.split('>')[0].strip(),
    url=url,
    license=license,
    keywords=["chemical kinetics", "Smoluchowski equation",
              "advection-diffusion-reaction"],
    packages=[pkg_name] + modules + tests,
    include_package_data=True,
    cmdclass=cmdclass_,
    ext_modules=ext_modules_,
    classifiers=classifiers,
    setup_requires=setup_requires,
    install_requires=install_requires,
    extras_require={'all': [
        'argh', 'pytest', 'scipy>=0.15', 'matplotlib', 'mpld3',
        'sym', 'sympy', 'pyodeint', 'pygslodeiv2', 'batemaneq',
        'sphinx', 'sphinx_rtd_theme', 'numpydoc'
    ]}
)

if __name__ == '__main__':
    try:
        if TAGGED_RELEASE:
            # Same commit should generate different sdist files
            # depending on tagged version (see RELEASE_VERSION)
            # this will ensure source distributions contain the correct version
            shutil.move(release_py_path, release_py_path+'__temp__')
            open(release_py_path, 'wt').write(
                "__version__ = '{}'\n".format(__version__))
        shutil.move(config_py_path, config_py_path+'__temp__')
        with open(config_py_path, 'wt') as fh:
            fh.write("env = {}\n".format(pprint.pformat(env)))
        setup(**setup_kwargs)
    finally:
        if TAGGED_RELEASE:
            shutil.move(release_py_path+'__temp__',
                        release_py_path)
        shutil.move(config_py_path+'__temp__',
                    config_py_path)
