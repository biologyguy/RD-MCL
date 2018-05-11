#! /bin/bash
# -*- coding: utf-8 -*-

export PATH="$HOME/miniconda/bin:$PATH"
FAILURE=0

# Disable py.test cacheprovider because it requires r/w access to the test directory
#### Pre-tests
cd /home/travis/build/biologyguy/RD-MCL/rdmcl/tests
TEST_SCRIPTS='test_fixtures.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Helper tests
TEST_SCRIPTS='test_helpers.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### MCMCMC tests
TEST_SCRIPTS='test_mcmcmc.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### merge_orthogroups tests
TEST_SCRIPTS='test_merge_orthogroups.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Worker tests
TEST_SCRIPTS='test_worker.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### RD-MCL tests
TEST_SCRIPTS='test_rdmcl.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Group by cluster tests
TEST_SCRIPTS='test_group_by_cluster.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Rename orthogroup tests
TEST_SCRIPTS='test_rename_orthogroup.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Reset workers tests
TEST_SCRIPTS='test_reset_workers.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Install tests
TEST_SCRIPTS='test_install.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Compare homolog groups tests
TEST_SCRIPTS='test_compare_homolog_groups.py '
py.test ${TEST_SCRIPTS} --cache-clear -p no:cacheprovider -p no:logging --cov --cov-append --cov-report= --cov-config ../.coveragerc --durations=10
if [ $? -ne 0 ]
then
    FAILURE=1
fi

#### Move coverage report file
mv .coverage /home/travis/build/biologyguy/RD-MCL/

#### Run Coveralls
cd /home/travis/build/biologyguy/RD-MCL
coveralls --config_file=/home/travis/build/biologyguy/RD-MCL/rdmcl/.coveragerc

exit ${FAILURE}