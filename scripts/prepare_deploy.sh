#!/bin/bash
touch doc/_build/html/.nojekyll
cp LICENSE doc/_build/html/
mkdir -p deploy/public_html/branches/"${CI_BRANCH}" deploy/script_queue
cp -r dist/* examples/ doc/_build/html/ deploy/public_html/branches/"${CI_BRANCH}"/
if bash -c '[[ "$CI_BRANCH" == "master" ]]'; then
    if [ -e benchmarks/ ]; then
        cat <<EOF>deploy/script_queue/run_benchmark.sh
source /etc/profile
cd ~/benchmarks/
git pull
asv run -k -e >asv-run.log
asv publish>asv-publish.log
EOF
        chmod +x deploy/script_queue/run_benchmark.sh
    fi
fi
