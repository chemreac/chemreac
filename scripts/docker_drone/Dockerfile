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
apt-get --quiet --assume-yes install git gcc g++ gfortran valgrind cmake wget curl liblapack-dev libsundials-serial-dev graphviz dot2tex && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install miniconda and install packages from binstar and PyPI
RUN wget --quiet http://repo.continuum.io/miniconda/Miniconda-${MINICONDA_VERSION}-Linux-x86_64.sh -O miniconda.sh && \
bash miniconda.sh -b -p $MINICONDA_PATH && \
rm miniconda.sh && \
export ORI_PATH=$PATH && \
export PATH="$MINICONDA_PATH/bin:$PATH" && \
conda config --set always_yes yes --set changeps1 no && \
conda config --add channels http://conda.binstar.org/bjodah && \
conda update --quiet conda && \
conda install --quiet conda-build binstar jinja2 numpy=1.9.1 scipy=0.14.1 matplotlib=1.4.2 cython=0.21.2 mako periodictable quantities pytest future pip sphinx numpydoc pycompilation pycodeexport sympy && \
pip install --quiet argh pytest-pep8 pytest-cov python-coveralls sphinx_rtd_theme mpld3 && \
rm -rf /tmp/* && \
mkdir -p $HOME/.config/matplotlib/ && echo "backend: Agg" > $HOME/.config/matplotlib/matplotlibrc && \
conda clean -t
