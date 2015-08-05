FROM ubuntu:trusty
MAINTAINER Björn Dahlgren <bjodah@DELETEMEgmail.com>
ENV DEBIAN_FRONTEND noninteractive
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV MINICONDA_VERSION 3.10.1
ENV MINICONDA_PATH /opt/miniconda
RUN locale-gen en_US.UTF-8

# Install packages from distribution repo
RUN apt-get update && \
apt-get --quiet --assume-yes install wget libsundials-serial graphviz dot2tex && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install miniconda and install packages from anaconda.org and PyPI
RUN wget --quiet http://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -O miniconda3.sh && \
bash miniconda3.sh -b -p $MINICONDA_PATH && \
rm miniconda3.sh && \
export ORI_PATH=$PATH && \
export PATH="$MINICONDA_PATH/bin:$PATH" && \
conda config --set always_yes yes --set changeps1 no && \
conda config --add channels http://conda.anaconda.org/bjodah && \
conda update --quiet conda && \
conda install --quiet numpy=1.9.2 scipy=0.16 matplotlib=1.4.3 periodictable quantities pytest future pip sympy chemreac && \
pip install --quiet argh pytest-pep8 pytest-cov python-coveralls mpld3 && \
rm -rf /tmp/* && \
mkdir -p /root/.config/matplotlib/ && echo "backend: Agg" > /root/.config/matplotlib/matplotlibrc && \
sed -i 's/Qt4Agg/Agg/g' /opt/miniconda/lib/python*/site-packages/matplotlib/mpl-data/matplotlibrc && \
sed -i 's/Qt4Agg/Agg/g' /opt/miniconda/pkgs/matplotlib-*/lib/*/site-packages/matplotlib/mpl-data/matplotlibrc && \
conda clean -t
