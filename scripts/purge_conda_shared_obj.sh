EXTMODULE=$1
CONDA_PATH=$2
SHARED_OBJ=$3

COUNTER=0
while true; do
    # Purge anaconda provided shared obj, e.g. libm.so problem:
    # see: https://groups.google.com/a/continuum.io/forum/#!msg/anaconda/4ysbSS2F8l4/2sHHioRKuMIJ
    CUR_SO=$(ldd $EXTMODULE | grep $SHARED_OBJ | awk '{print $3}')
    if [[ "$CUR_SO" == "${CONDA_PATH}"* ]]; then
        mv $CUR_SO _$(basename ${CUR_SO}).${COUNTER}
        COUNTER=$((COUNTER + 1))
    else
        echo $CUR_SO
        break
    fi
done
