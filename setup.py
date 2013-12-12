#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from distutils.core import setup
from pycompilation.dist import clever_build_ext, CleverExtension


setup(
    name='chemreac',
    version='0.0.1',
    description='Python extension for reaction diffusion.',
    author='Björn Dahlgren',
    author_email='bjodah@DELETEMEgmail.com',
    #url='https://github.com/bjodah/###',
    packages=['chemreac'],
    cmdclass = {'build_ext': clever_build_ext},
    ext_modules = [
        CleverExtension(
            "chemreac.cpp_chem_wrapper",
            sources=[
                'chemreac/cpp_chem.cpp',
                'chemreac/cpp_chem_wrapper.pyx',
            ],
            template_regexps=[
                (r'(\w+)_template(\w+)', r'\1\2', {}),
            ]
    ]
)
