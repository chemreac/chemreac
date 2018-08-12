#!/bin/bash -x
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
    echo ${CI_BRANCH} | tail -c +2 > __conda_version__.txt
fi

python3 -m pip install argh finitediff block_diag_ilu pycvodes

set -e

(cd tests-native; make -B CONTEXT=valgrind EXTRA_COMPILE_ARGS='-D_GLIBCXX_DEBUG' test)
(cd tests-native; make -B CXX=clang++-6.0 CC=clang-6.0 OPTIMIZE=1 WITH_OPENMP=0 EXTRA_COMPILE_ARGS='-fsanitize=address -DNDEBUG' test)

python3 -m pip install -e .[all]
git clean -xfd
CFLAGS="-D_GLIBCXX_DEBUG" python3 setup.py develop
bash -c "ulimit -v 2048000; ./scripts/run_tests.sh"
git clean -xfd
CC=clang-6.0 \
  CXX=clang++-6.0 \
  CFLAGS="-fsanitize=address -UNDEBUG" \
  python3 setup.py build_ext -i
ASAN_OPTIONS=detect_leaks=0 LD_PRELOAD=/usr/lib/llvm-6.0/lib/clang/6.0.1/lib/linux/libclang_rt.asan-x86_64.so ./scripts/run_tests.sh ${@:2}
git clean -xfd

python3 setup.py sdist
cp dist/${PKG_NAME}-*.tar.gz /tmp
(cd /; python3 -m pip install --ignore-installed /tmp/${PKG_NAME}-*.tar.gz; python3 -c "import $PKG_NAME")
