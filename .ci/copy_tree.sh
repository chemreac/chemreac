if [ "$#" -ne 2 ]; then
    2>&1 echo "Two arguments needed (destination folder, sundials base)"
    exit 1
fi
cp -ra ci_cache/pyusrb /opt
export PYTHONUSERBASE=/opt/pyusrb

git-archive-all --prefix='' /tmp/HEAD.zip

set -u
DESTDIR=$1
SUNDBASE=$2
set +u
export CFLAGS="-isystem $SUNDBASE/include $CFLAGS"
export LDFLAGS="-Wl,--disable-new-dtags -Wl,-rpath,$SUNDBASE/lib -L$SUNDBASE/lib $LDFLAGS"

mkdir -p "$DESTDIR"
cd "$DESTDIR"
unzip /tmp/HEAD.zip

${PYTHON:-python3} -c "import pycvodes as pcv; print(pcv.__version__)" || ( >&2 echo "failed to install pycvodes"; exit 1 )
${PYTHON:-python3} -c "import block_diag_ilu as bdi; print(bdi.__version__)" || ( >&2 echo "failed to install block_diag_ilu"; exit 1 )
