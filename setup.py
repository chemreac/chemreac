#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

from distutils.core import setup

name_ = 'chemreac'
version_ = '0.0.4'

DEBUG=True
USE_OPENMP = os.environ.get('USE_OPENMP', False)

if '--help'in sys.argv[1:] or sys.argv[1] in (
        '--help-commands', 'egg_info', 'clean', '--version'):
    cmdclass_ = {}
    ext_modules_ = []
else:
    from pycompilation.dist import clever_build_ext, CleverExtension
    cmdclass_ = {'build_ext': clever_build_ext}
    ext_modules_ = [
        CleverExtension(
            "chemreac._chemreac",
            sources=[
                'src/chemreac_template.cpp',
                'chemreac/_chemreac.pyx',
            ],
            template_regexps=[
                (r'^(\w+)_template.(\w+)$', r'\1.\2', {'USE_OPENMP': USE_OPENMP}),
            ],
            pycompilation_compile_kwargs={
                'per_file_kwargs': {
                    'src/chemreac.cpp': {
                        'std': 'c++11',
                        'options': ['pic', 'warn', 'fast'] +\
                        (['openmp'] if USE_OPENMP else []),
                        'defmacros': ['restrict=__restrict__']+\
                        (['DEBUG'] if DEBUG else []),
                    },
                },
                'defmacros': ['restrict=__restrict__'],
            },
            pycompilation_link_kwargs={
                'options': (['openmp'] if USE_OPENMP else []),
            },
            include_dirs=['src/'],
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
    packages=[name_],
    cmdclass = cmdclass_,
    ext_modules = ext_modules_,
)
