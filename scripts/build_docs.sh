#!/bin/bash -e
cd docs/
PYTHONPATH=`pwd`/..:`pwd`/../examples make html >_build.log
cp _build.log _build/html/
