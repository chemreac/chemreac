#!/bin/bash
SUNDIALS_URLS=("http://computation.llnl.gov/casc/sundials/download/code/sundials-2.5.0.tar.gz" "http://pkgs.fedoraproject.org/repo/pkgs/sundials/sundials-2.5.0.tar.gz/aba8b56eec600de3109cfb967aa3ba0f/sundials-2.5.0.tar.gz" "http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sundials-2.5.0.tar.gz")

SUNDIALS_FNAME="sundials-2.5.0.tar.gz"
SUNDIALS_MD5="aba8b56eec600de3109cfb967aa3ba0f"

for URL in "${SUNDIALS_URLS[@]}"; do
    echo "Downloading ${URL}..."
    wget -q $URL
    if [ $? -eq 0 ]; then
        DOWNLOAD_MD5=$(md5sum ${SUNDIALS_FNAME} | cut -d ' ' -f 1)
        echo "md5: ${DOWNLOAD_MD5}"
        if [[ "$DOWNLOAD_MD5" == "$SUNDIALS_MD5" ]]; then
            tar xzf $SUNDIALS_FNAME
            mkdir sundials_build
            cd sundials_build
            cmake -DBUILD_SHARED_LIBS:BOOL="1" -DCMAKE_BUILD_TYPE:STRING="Debug" -DLAPACK_ENABLE:BOOL=1 ../sundials-*/
            make
            sudo make install
            exit 0
        fi
    fi    
done
exit 1
