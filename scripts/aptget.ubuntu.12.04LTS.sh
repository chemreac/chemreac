#!/bin/bash
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get --quiet update
sudo apt-get --quiet --assume-yes install gfortran-4.8 g++-4.8 gcc-4.8 liblapack-dev cmake graphviz dot2tex valgrind
# These are handled by miniconda:
#python-numpy python-scipy cython

if [[ $? != 0 ]]; then
    echo "apt-get install failed."
    exit 1
fi
sudo update-alternatives --remove-all gcc 
sudo update-alternatives --remove-all g++
sudo update-alternatives --remove-all gfortran
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 20
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 20
sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-4.8 20
sudo update-alternatives --config gcc
sudo update-alternatives --config g++
sudo update-alternatives --config gfortran
gcc --version
g++ --version
gfortran --version
