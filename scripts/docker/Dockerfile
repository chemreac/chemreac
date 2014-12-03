FROM ubuntu:trusty
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get --quiet --assume-yes install gcc g++ gfortran liblapack-dev libsundials-serial-dev graphviz dot2tex python-pip python-numpy python-scipy python-argh python-mako python-quantities python-matplotlib wget python-dev
RUN pip install pytest pytest-pep8 periodictable future mpld3 sympy pycodeexport pycompilation pytest-cov
ADD ./run-tests.sh /tmp/run-tests.sh
ENTRYPOINT ["/tmp/run-tests.sh"]