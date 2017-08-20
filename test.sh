#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
FAILURE=0

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
cd /home/travis/build/biologyguy/RD-MCL/rdmcl/tests
TEST_SCRIPTS='test_fixtures.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Helper tests
TEST_SCRIPTS='test_helpers.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Worker tests
TEST_SCRIPTS='test_worker.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### RD-MCL tests
TEST_SCRIPTS='test_rdmcl.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi
mv .coverage /home/travis/build/biologyguy/RD-MCL/


#### Run Coveralls
cd /home/travis/build/biologyguy/RD-MCL
coveralls

exit ${FAILURE}