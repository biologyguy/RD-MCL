#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "test_fixtures.py"
TEST_SCRIPTS=${DIR}'/tests/test_fixtures.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/tests/__init__.py --cov-report html -n 3 -p no:cacheprovider --durations=10 $@

echo "test_rdmcl.py"
TEST_SCRIPTS=${DIR}'/tests/test_rdmcl.py '
py.test ${TEST_SCRIPTS} --cov ${DIR}/rdmcl.py --cov-report html -n 8 -p no:cacheprovider --durations=10 $@
