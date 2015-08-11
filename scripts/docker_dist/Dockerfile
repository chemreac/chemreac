FROM bjodah/dockerfile-bjodah-scicomp
MAINTAINER Björn Dahlgren <bjodah@DELETEMEgmail.com>
# FROM ubuntu:trusty
# ENV DEBIAN_FRONTEND noninteractive
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8
RUN locale-gen en_US.UTF-8
RUN wget --quiet "http://hera.physchem.kth.se/~chemreac/chemreac-0.3.0.tar.gz" -O /resources/ && \
export CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH:/usr/include/python2.7" && \
/scripts/install_py_tar_gz.sh chemreac-0.3.0
