#!/bin/bash
BUILD_TYPE=${1:-Debug}
SUNDIALS_URLS=(\
"http://hera.physchem.kth.se/~bjorn/sundials-2.5.0.tar.gz" \
"http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sundials-2.5.0.tar.gz" \
"http://pkgs.fedoraproject.org/repo/pkgs/sundials/sundials-2.5.0.tar.gz/aba8b56eec600de3109cfb967aa3ba0f/sundials-2.5.0.tar.gz" \
"http://computation.llnl.gov/casc/sundials/download/code/sundials-2.5.0.tar.gz" \
)

SUNDIALS_FNAME="sundials-2.5.0.tar.gz"
SUNDIALS_MD5="aba8b56eec600de3109cfb967aa3ba0f"

for URL in "${SUNDIALS_URLS[@]}"; do
    echo "Downloading ${URL}..."
    ( wget --quiet --tries=2 --timeout=45 $URL & sleep 60; kill %1 )
    if [ $? -eq 0 ]; then
        DOWNLOAD_MD5=$(md5sum ${SUNDIALS_FNAME} | cut -d ' ' -f 1)
        echo "md5: ${DOWNLOAD_MD5}"
        if [[ "$DOWNLOAD_MD5" != "$SUNDIALS_MD5" ]]; then
            rm ${SUNDIALS_FNAME}
        else
            tar xzf $SUNDIALS_FNAME
            mkdir sundials_build
            cd sundials_build
            cmake -DBUILD_SHARED_LIBS:BOOL="1" -DCMAKE_BUILD_TYPE:STRING="$1" -DLAPACK_ENABLE:BOOL=1 ../sundials-*/
            if [[ $? != 0 ]]; then
                echo "cmake of sundials failed."
                exit 1
            fi
            make
            if [[ $? != 0 ]]; then
                echo "make of sundials failed."
                exit 1
            fi
            sudo make install
            if [[ $? != 0 ]]; then
                echo "make install of sundials failed."
                exit 1
            fi
            exit 0
        fi
    fi    
done
exit 1
