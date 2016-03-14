#!/bin/bash -xe

ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
( cd $ABS_REPO_PATH/docs/; PYTHONPATH=$ABS_REPO_PATH:$ABS_REPO_PATH/examples make $@ html | tee _build.log; cp _build.log _build/html/ )
