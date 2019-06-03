if [ "$#" -ne 3 ]; then
    2>&1 echo "Three arguments needed (destination folder, symengine folder, sundials base)"
    exit 1
fi
mkdir -p "$1"
cp -ra external scripts kinetgine tests dist py-ext-mod util *.* .ci "$1"
cp -ra ci-cache/pyusrb /opt
export SUNDBASE=$2
export CPATH=$SUNDBASE/include:$CPATH
export LIBRARY_PATH=$SUNDBASE/lib
export LD_LIBRARY_PATH=$SUNDBASE/lib
export CMAKE_PREFIX_PATH=$SUNDBASE:$3:$CMAKE_PREFIX_PATH
cd $1
${PYTHON:-python3} -c "import pycvodes"
