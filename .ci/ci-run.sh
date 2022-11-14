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
( cd dist; python3 -c "import chemreac" )

(cd tests-native; make -B CONTEXT=valgrind EXTRA_COMPILE_ARGS='-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC' test)
(cd tests-native; make -B CXX=clang++-12 CC=clang-12 OPTIMIZE=1 WITH_OPENMP=0 EXTRA_COMPILE_ARGS='-fsanitize=address -DNDEBUG' test)
(cd tests-native; make -B CXX=clang++-12 CC=clang-12 OPTIMIZE=1 WITH_OPENMP=0 EXTRA_COMPILE_ARGS='-fsanitize=undefined' test)

CFLAGS="-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC $CFLAGS" python3 setup.py build_ext -i
bash -c "ulimit -v 3072000; ./scripts/run_tests.sh ${@:2}"

for san in address undefined; do
    if [[ $san == "address" ]]; then
        PY_LD_PRELOAD=$(clang++-12 --print-file-name=libclang_rt.asan-$( uname -m).so)
        SAN_OPTS="ASAN_OPTIONS=abort_on_error=1:detect_leaks=0"
    elif [[ $san == "undefined" ]]; then
        PY_LD_PRELOAD=$(clang++-12 --print-file-name=libclang_rt.ubsan_standalone-$(uname -m).so)
        SAN_OPTS="UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1"
    else
        >&2 echo "bug in script"
        exit 1
    fi
    rm -rf build/
    CC=clang-12 \
      CXX=clang++-12 \
      CFLAGS="-fsanitize=$san -UNDEBUG $CFLAGS" \
      python3 setup.py build_ext -i
    LD_PRELOAD=$PY_LD_PRELOAD \
              PYTHONMALLOC=malloc \
              env $SAN_OPTS \
              ./scripts/run_tests.sh -s -v 
done

python3 -m pip uninstall -y ${PKG_NAME}

rm -rf build/

python3 setup.py sdist
cp dist/${PKG_NAME}-*.tar.gz /tmp
(mkdir /tmp/sdist_tar_gz; cd /tmp/sdist_tar_gz; tar xf ../${PKG_NAME}-*.tar.gz; cd ${PKG_NAME}-*/; python3 setup.py build_ext -i; PYTHONPATH=$(pwd); python3 -c "import $PKG_NAME")

# Make sure repo is pip installable from git-archive zip
(mkdir /tmp/archive_zip; cd /tmp/archive_zip; unzip -q ../HEAD.zip; python3 setup.py build_ext -i; PYTHONPATH=$(pwd); python3 -c "import $PKG_NAME")