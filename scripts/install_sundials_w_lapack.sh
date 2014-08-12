#!/bin/bash
SUNDIALS_URLS=("http://pkgs.fedoraproject.org/repo/pkgs/sundials/sundials-2.5.0.tar.gz/aba8b56eec600de3109cfb967aa3ba0f/sundials-2.5.0.tar.gz" "http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sundials-2.5.0.tar.gz")

SUNDIALS_FNAME="sundials-2.5.0.tar.gz"
SUNDIALS_MD5="4459e74e19820e7a5c6fe2dc1b0c2912"

for URL in "${SUNDIALS_URLS[@]}"; do
    wget $URL
    if [ $? -eq 0 ]; then
        DOWNLOAD_MD5 = $(md5sum ${SUNDIALS_FNAME} | cut -d ' ' -f 1)
        if [ "$DOWNLOAD_MD5" == "$SUNDIALS_MD5" ]; then
            tar xvzf $SUNDIALS_FNAME
            mkdir sundials_build
            cd sundials_build
            cmake --enable-shared --enable-lapack ../sundials-*/
            make
            sudo make install
            exit 0
        fi
    fi    
done
exit 1
