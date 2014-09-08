PY_VERSION=$1
INSTALL_PATH=$2

if [[ "$PY_VERSION" == "2.7" ]]; then
    # This saves us some downloading for this version
    wget --quiet http://repo.continuum.io/miniconda/Miniconda-3.6.0-Linux-x86_64.sh -O miniconda.sh;
else
    wget --quiet http://repo.continuum.io/miniconda/Miniconda3-3.6.0-Linux-x86_64.sh -O miniconda.sh;
fi
if [[ $? != 0 ]]; then
    echo "Failed to get Miniconda."
    exit 1
fi

chmod +x miniconda.sh
./miniconda.sh -b -p $INSTALL_PATH
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update --quiet conda
conda info -a
