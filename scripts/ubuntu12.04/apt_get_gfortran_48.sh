#!/bin/bash
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get --quiet update
sudo apt-get --quiet --assume-yes install gfortran-4.8
sudo update-alternatives --remove-all gfortran
sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-4.8 20
sudo update-alternatives --config gfortran
GFORTRAN_VERSION=$(gfortran --version | grep -i fortran | head -1 | sed 's/^.* //g')
if [ "${GFORTRAN_VERSION:0:3}" != "4.8"  ]; then
    echo "Failed to install gfortran-4.8"
    exit 1
fi
