#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

from distutils.core import setup, Command

name_ = 'chemreac'
version_ = '0.2.1'

DEBUG = True if os.environ.get('USE_DEBUG', False) else False
USE_OPENMP = True if os.environ.get('USE_OPENMP', False) else False
LLAPACK = os.environ.get('LLAPACK', 'lapack')

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
on_drone = os.environ.get('DRONE', 'false') == 'true'
on_travis = os.environ.get('TRAVIS', 'flse') == 'true'

if on_drone or on_travis:
    # 'fast' implies march=native which fails on current version of docker.
    options = ['pic', 'warn']
else:
    options = ['pic', 'warn', 'fast']

# Make `python setup.py test` work without depending on py.test being installed
# https://pytest.org/latest/goodpractises.html


class PyTest(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        import sys
        # py.test --genscript=runtests.py
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

cmdclass_ = {'test': PyTest}

if on_rtd or '--help' in sys.argv[1:] or sys.argv[1] in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    # Enbale pip to probe setup.py before all requirements are installed
    ext_modules_ = []
else:
    import pickle
    from pycodeexport.dist import pce_build_ext, PCEExtension
    import numpy as np
    cmdclass_['build_ext'] = pce_build_ext
    subsd = {'USE_OPENMP': USE_OPENMP}
    sources = [
        'src/chemreac_template.cpp',
        'src/finitediff/finitediff/fornberg.f90',
        'src/finitediff/finitediff/c_fornberg.f90',
        'chemreac/_chemreac.pyx',
    ]

    ext_modules_ = [
        PCEExtension(
            "chemreac._chemreac",
            sources=sources,
            template_regexps=[
                (r'^(\w+)_template.(\w+)$', r'\1.\2', subsd),
            ],
            pycompilation_compile_kwargs={
                'per_file_kwargs': {
                    'src/chemreac.cpp': {
                        'std': 'c++0x',
                        # 'fast' doesn't work on drone.io
                        'options': options +
                        (['openmp'] if USE_OPENMP else []),
                        'defmacros': ['DEBUG'] +
                        (['DEBUG'] if DEBUG else []),
                    },
                    'src/chemreac_sundials.cpp': {
                        'std': 'c++0x',
                        'options': options
                    },
                },
                'options': options,
            },
            pycompilation_link_kwargs={
                'options': (['openmp'] if USE_OPENMP else []),
                'std': 'c++0x',
                'libs': ['sundials_cvode', LLAPACK, 'sundials_nvecserial'],
            },
            include_dirs=['src/', 'src/finitediff/finitediff/',
                          np.get_include()],
            logger=True,
        )
    ]

setup(
    name=name_,
    version=version_,
    description='Python extension for reaction diffusion.',
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://bitbucket.org/bjodah/'+name_,
    packages=[name_, name_+'.util'],
    cmdclass=cmdclass_,
    ext_modules=ext_modules_,
)
