#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "test_rdmcl.py"
TEST_SCRIPTS=${DIR}'/tests/test_rdmcl.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/rdmcl.py --cov-report html -n 1 -p no:cacheprovider --durations=10 $@
