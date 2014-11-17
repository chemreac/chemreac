export CONDA_PATH=$(./scripts/get_conda_path.sh)
#mv ${CONDA_PATH}/lib/libm.so libm_.so_
#./scripts/purge_conda_shared_obj.sh $(find . -name _chemreac.so | head -1) $CONDA_PATH libm
$PYTHON setup.py install
