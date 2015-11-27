export LIBRARY_PATH="$PREFIX/lib:$LIBRARY_PATH"
export INCLUDE_PATH="$PREFIX/include:$INCLUDE_PATH"
$PYTHON setup.py build
$PYTHON setup.py install
