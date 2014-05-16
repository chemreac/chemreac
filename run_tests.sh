#!/bin/bash
py.test
echo `git rev-parse HEAD` `hostname` `grep Best src/tests/test_chemreac.out | cut -d: -f2` >>performance_log.txt
