#!/bin/bash
cd docs/
PYTHONPATH=`pwd`/.. make html >_build.log
cp _build.log _build/html/
