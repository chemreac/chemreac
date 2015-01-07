Trusty for CI testing of chemreac using drone
=============================================

Comes with miniconda (inspect version $MINICONDA_VERSION) installed
in $MINICONDA_PATH, $PATH also has $MINICONDA_PATH/bin in it.

This image essentially provides caching for the CI so that
all dependencies don't need to be downloaded for each test run.

See Dockerfile for installed packages.
