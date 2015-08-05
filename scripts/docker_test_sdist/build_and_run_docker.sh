#!/bin/bash
docker build -t chemreac-sdist .
docker run -t chemreac-sdist ${1:-chemreac-0.3.rc7}
