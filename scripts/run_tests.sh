#!/bin/bash
py.test --slow --veryslow --pep8 --doctest-modules --cov chemreac --cov-report html --ignore setup.py
