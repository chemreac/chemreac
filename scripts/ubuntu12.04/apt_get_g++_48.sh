#!/bin/bash
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get --quiet update
sudo apt-get --quiet --assume-yes install g++-4.8
sudo update-alternatives --remove-all g++
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 20
sudo update-alternatives --config g++
GPP_VERSION=$(g++ --version | grep ^g++ | sed 's/^.* //g')
if [[ "${GPP_VERSION:0:3}" -ne "4.8"  ]]; then
    echo "Failed to install g++-4.8"
    exit 1
fi
