export LIBRARY_PATH="$PREFIX/lib:$LIBRARY_PATH"
export CPLUS_INCLUDE_PATH="$PREFIX/include:$INCLUDE_PATH"
$PYTHON setup.py build
$PYTHON setup.py install
