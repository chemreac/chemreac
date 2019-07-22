#!/bin/bash -xe

ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
export MPLBACKEND=Agg
( cd $ABS_REPO_PATH/doc/; PYTHONPATH=$ABS_REPO_PATH:$ABS_REPO_PATH/examples:$PYTHONPATH make $@ html | tee _build.log; cp _build.log _build/html/ )
