#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from distutils.core import setup
from pycompilation.dist import clever_build_ext, CleverExtension

name_ = 'chemreac'
version_ = '0.0.1'

USE_OPENMP = os.environ.get('USE_OPENMP', False)

setup(
    name=name_,
    version=version_,
    description='Python extension for reaction diffusion.',
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    url='https://bitbucket.org/bjodah/'+name_,
    packages=[name_],
    cmdclass = {'build_ext': clever_build_ext},
    ext_modules = [
        CleverExtension(
            "chemreac.cpp_chem_wrapper",
            sources=[
                'chemreac/cpp_chem_template.cpp',
                'chemreac/cpp_chem_wrapper.pyx',
            ],
            template_regexps=[
                (r'^(\w+)_template.(\w+)$', r'\1.\2', {'USE_OPENMP': USE_OPENMP}),
            ],
            pycompilation_compile_kwargs={
                'per_file_kwargs': {
                    'chemreac/cpp_chem.cpp': {
                        'std': 'c++11',
                        'options': ['pic', 'warn', 'fast'] +\
                        (['openmp'] if USE_OPENMP else [])
                    },
                },
            },
            pycompilation_link_kwargs={
                'options': (['openmp'] if USE_OPENMP else []),
            },
            include_dirs=['chemreac/'],
            define_macros=['DEBUG', 'restrict=__restrict__'],
            logger=True,
        )
    ]
)
