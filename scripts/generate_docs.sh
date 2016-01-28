#!/bin/bash -xe

ABS_REPO_PATH=$(unset CDPATH && cd "$(dirname "$0")/.." && echo $PWD)
( cd docs/; PYTHONPATH=$ABS_REPO_PATH:$ABS_REPO_PATH/examples make $@ html >_build.log; cp _build.log _build/html/ )
