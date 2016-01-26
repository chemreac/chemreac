#!/bin/bash -xe

ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
export PYTHONPATH=$absolute_repo_path:$absolute_repo_path/examples:$PYTHONPATH
( cd docs/; PYTHONPATH=$ABS_REPO_PATH make $@ html >_build.log; cp _build.log _build/html/ )
