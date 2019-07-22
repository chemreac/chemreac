#!/bin/bash -x
PKG_NAME=${1:-${CI_REPO##*/}}
if [[ "$CI_BRANCH" =~ ^v[0-9]+.[0-9]?* ]]; then
    eval export ${1^^}_RELEASE_VERSION=\$CI_BRANCH
fi

python3 -m pip install --user -e .[all]
python3 setup.py build_ext -i
mkdir -p dist
cp -r chemreac dist/.
ls dist/chemreac
rm -r build
set -e

(cd tests-native; make -B CONTEXT=valgrind EXTRA_COMPILE_ARGS='-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC' test)
(cd tests-native; make -B CXX=clang++-8 CC=clang-8 OPTIMIZE=1 WITH_OPENMP=0 EXTRA_COMPILE_ARGS='-fsanitize=address -DNDEBUG' test)

CFLAGS="-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" python3 setup.py build_ext -i
bash -c "ulimit -v 2048000; ./scripts/run_tests.sh"

rm -rf build/
CC=clang-8 \
  CXX=clang++-8 \
  CFLAGS="-fsanitize=address -UNDEBUG" \
  python3 setup.py build_ext -i
PYTHONMALLOC=malloc ASAN_OPTIONS=detect_leaks=0 LD_PRELOAD=/usr/lib/llvm-8/lib/clang/8.0.1/lib/linux/libclang_rt.asan-x86_64.so ./scripts/run_tests.sh "${@:2}"

python3 -m pip uninstall -y chemreac

rm -rf build/

python3 setup.py sdist
cp dist/${PKG_NAME}-*.tar.gz /tmp
(mkdir /tmp/sdist_tar_gz; cd /tmp/sdist_tar_gz; tar xf ../${PKG_NAME}-*.tar.gz; cd ${PKG_NAME}-*/; python3 setup.py build_ext -i; PYTHONPATH=$(pwd); python3 -c "import $PKG_NAME")

# Make sure repo is pip installable from git-archive zip
(mkdir /tmp/archive_zip; cd /tmp/archive_zip; unzip ../HEAD.zip; python3 setup.py build_ext -i; PYTHONPATH=$(pwd); python3 -c "import $PKG_NAME")
