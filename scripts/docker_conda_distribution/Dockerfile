FROM bjodah/dockerfile-bjodah-scicomp
MAINTAINER Björn Dahlgren <bjodah@DELETEMEgmail.com>
# FROM ubuntu:trusty
# ENV DEBIAN_FRONTEND noninteractive
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
ENV CHEMREAC_VERSION 0.3.0
RUN locale-gen en_US.UTF-8

# Install packages from distribution repo
RUN apt-get update && \
apt-get --quiet --assume-yes install wget libsundials-serial graphviz dot2tex && \
apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install miniconda and install packages from anaconda.org and PyPI
RUN wget --quiet "http://hera.physchem.kth.se/~chemreac/chemreac-${CHEMREAC_VERSION}.tar.gz" -O chemreac.tar.gz && \
tar xvzf chemreac.tar.gz && \
rm chemreac.tar.gz && \
cd chemreac*/ && \
python setup.py install && \
cd - && \
rm -r chemreac*/ && \
py.test -m --pyargs chemreac