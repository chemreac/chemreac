if [ "$#" -ne 2 ]; then
    2>&1 echo "Two arguments needed (destination folder, sundials base)"
    exit 1
fi
mkdir -p "$1"
cp -ra external scripts chemreac tests-native examples *.* .ci "$1"
cp -ra ci-cache/pyusrb /opt
export PYTHONUSERBASE=/opt/pyusrb

export SUNDBASE=$2
export CPATH=$SUNDBASE/include:$CPATH
export LIBRARY_PATH=$SUNDBASE/lib
export LD_LIBRARY_PATH=$SUNDBASE/lib
export CMAKE_PREFIX_PATH=$SUNDBASE:$CMAKE_PREFIX_PATH
cd $1
${PYTHON:-python3} -c "import pycvodes"
