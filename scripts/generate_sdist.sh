#!/bin/bash -e
cd $(dirname $0)/..
python setup.py sdist
scp ./dist/$(ls -rt ./dist | tail -1) chemreac@hera.physchem.kth.se:~/public_html/
ssh chemreac@hera.physchem.kth.se '~/public_html/generate_sdist_index.sh'
