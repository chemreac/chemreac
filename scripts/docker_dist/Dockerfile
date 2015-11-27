FROM bjodah/bjodahimg:v1.0
MAINTAINER Björn Dahlgren <bjodah@DELETEMEgmail.com>
RUN wget --quiet "http://hera.physchem.kth.se/~chemreac/chemreac-0.3.0.tar.gz" -O /tmp/ && \
export CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH:/usr/include/python2.7" && \
cd /tmp && tar xzf chemreac-0.3.0.tar.gz && cd chemreac && python setup.py install
