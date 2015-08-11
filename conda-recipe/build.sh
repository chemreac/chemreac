export LIBRARY_PATH="$PREFIX/lib"
export INCLUDE_PATH="$PREFIX/include"
echo "LIBRARY_PATH=${LIBRARY_PATH}"
ls ${LIBRARY_PATH}
$PYTHON setup.py install
