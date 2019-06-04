if [ "$#" -ne 2 ]; then
    2>&1 echo "Two arguments needed (destination folder, sundials base)"
    exit 1
fi
mkdir -p "$1"
git archive HEAD | tar x -C "$1"
cp -ra ci_cache/pyusrb /opt
export PYTHONUSERBASE=/opt/pyusrb

export SUNDBASE=$2
export CPATH=$SUNDBASE/include:$CPATH
export LIBRARY_PATH=$SUNDBASE/lib
export LD_LIBRARY_PATH=$SUNDBASE/lib
export CMAKE_PREFIX_PATH=$SUNDBASE:$CMAKE_PREFIX_PATH
python3 -m pip install --user "pycvodes<0.11.0"
cd $1
${PYTHON:-python3} -c "import pycvodes as pcv; print(pcv.__version__)" || ( >&2 echo "failed to install pycvodes"; exit 1 )
${PYTHON:-python3} -c "import block_diag_ilu as bdi; print(bdi.__version__)" || ( >&2 echo "failed to install block_diag_ilu"; exit 1 )
