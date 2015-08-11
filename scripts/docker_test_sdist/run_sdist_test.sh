#!/bin/bash
SLUG=${1:-chemreac-0.3.0}
TMPD=$(mktemp -d)
SCRIPTNAME=my_script.sh
FULLPATH=$TMPD/$SCRIPTNAME
cat <<EOF >$FULLPATH
#!/bin/bash
wget "http://hera.physchem.kth.se/~chemreac/$SLUG.tar.gz"
tar xvzf "$SLUG.tar.gz"
pushd "$SLUG"
CPLUS_INCLUDE_PATH=/usr/include/python2.7 ./scripts/run_tests.sh
popd
EOF
chmod +x $FULLPATH
trap "rm -r $TMPD; docker rm chemreac_sdist_test"
docker run --name chemreac_sdist_test -v $TMPD:/input:ro -w /tmp -i bjodah/bjodah-scicomp /input/$SCRIPTNAME
docker wait chemreac_sdist_test
