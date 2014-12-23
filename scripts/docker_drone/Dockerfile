FROM ubuntu:trusty
ENV DEBIAN_FRONTEND noninteractive
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV MINICONDA_VERSION 3.7.3
ENV MINICONDA_PATH /opt/miniconda
RUN locale-gen en_US.UTF-8

# Install packages from distribution repo
RUN apt-get update && \
apt-get --quiet --assume-yes install gcc g++ gfortran valgrind cmake wget curl liblapack-dev libsundials-serial-dev graphviz dot2tex python-pip python-numpy python-scipy python-argh python-mako python-quantities python-matplotlib python-dev && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install python packages from PyPI
RUN pip install pytest pytest-pep8 periodictable future mpld3 sympy pycodeexport pycompilation pytest-cov sphinx && \
rm -rf /tmp/* && \
mkdir -p $HOME/.config/matplotlib/ && echo "backend: Agg" > $HOME/.config/matplotlib/matplotlibrc

# Install miniconda
RUN wget --quiet http://repo.continuum.io/miniconda/Miniconda-${MINICONDA_VERSION}-Linux-x86_64.sh -O miniconda.sh && \
bash miniconda.sh -b -p $MINICONDA_PATH && \
rm miniconda.sh && \
export ORI_PATH=$PATH && \
export PATH="$MINICONDA_PATH/bin:$PATH" && \
conda config --set always_yes yes --set changeps1 no && \
PATH="$MINICONDA_PATH/bin:$PATH" conda update --quiet conda && \
PATH="$MINICONDA_PATH/bin:$PATH" conda install --quiet conda-build binstar jinja2 cython && \
PATH="$MINICONDA_PATH/bin:$PATH" conda clean -t && \
export PATH=$ORI_PATH

# Oops, forgot about git
RUN apt-get update && apt-get install --quiet --assume-yes git