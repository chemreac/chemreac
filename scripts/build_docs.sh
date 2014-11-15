#!/bin/bash -e

# Extract absolute path of script, from:
# http://unix.stackexchange.com/a/9546
# Note: we are assuming this script is inside a subdirectory of the repo root
absolute_repo_path_x="$(readlink -fn -- "$(dirname $0)/.."; echo x)"
absolute_repo_path="${absolute_repo_path_x%x}"
export PYTHONPATH=$absolute_repo_path:$absolute_repo_path/examples:$PYTHONPATH

cd docs/
make $@ html >_build.log
cp _build.log _build/html/
