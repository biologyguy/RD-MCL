#!/usr/bin/env bash

find . -name "__pycache__" -type d | xargs rm -r || echo "No pycache detected"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "test_fixtures.py"
TEST_SCRIPTS=${DIR}'/tests/test_fixtures.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.tests.__init__ --cov-report html --cov-config .coveragerc -n 4 -p no:cacheprovider -p no:logging --durations=10 $@

echo "test_helpers.py"
TEST_SCRIPTS=${DIR}'/tests/test_helpers.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.helpers --cov-report html --cov-config .coveragerc -n 8 -p no:cacheprovider -p no:logging --durations=10 $@

echo "test_mcmcmc.py"
TEST_SCRIPTS=${DIR}'/tests/test_mcmcmc.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.mcmcmc --cov-report html --cov-config .coveragerc -n 8 -p no:cacheprovider -p no:logging --durations=10 $@

echo "test_merge_orthogroups.py"
TEST_SCRIPTS=${DIR}'/tests/test_merge_orthogroups.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.merge_orthogroups --cov-report html -n 8 --cov-config .coveragerc -p no:cacheprovider -p no:logging --durations=10 -m "not slow" $@

echo "test_worker.py"
TEST_SCRIPTS=${DIR}'/tests/test_worker.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.launch_worker --cov-report html --cov-config .coveragerc -n 8 -p no:cacheprovider -p no:logging --durations=10 $@

echo "test_group_by_cluster.py"
TEST_SCRIPTS=${DIR}'/tests/test_group_by_cluster.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.group_by_cluster --cov-report html --cov-config .coveragerc -n 8 -p no:cacheprovider -p no:logging --durations=10 $@

echo "test_rdmcl.py"
TEST_SCRIPTS=${DIR}'/tests/test_rdmcl.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.rdmcl --cov-report html -n 8 --cov-config .coveragerc -p no:cacheprovider -p no:logging --durations=10 -m "not slow" $@

echo "test_rename_orthogroup.py"
TEST_SCRIPTS=${DIR}'tests/test_rename_orthogroup.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.rename_orthogroup --cov-report html -n 8 --cov-config .coveragerc -p no:cacheprovider -p no:logging --durations=10  $@

echo "test_reset_workers.py"
TEST_SCRIPTS=${DIR}'tests/test_reset_workers.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.reset_workers --cov-report html -n 8 --cov-config .coveragerc -p no:cacheprovider -p no:logging --durations=10 -m $@

echo "test_install.py"
TEST_SCRIPTS=${DIR}'tests/test_install.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.install --cov-report html -n 8 --cov-config .coveragerc -p no:cacheprovider -p no:logging --durations=10 -m $@

echo "test_compare_homolog_groups.py"
TEST_SCRIPTS=${DIR}'tests/test_compare_homolog_groups.py '
py.test ${TEST_SCRIPTS} --cov rdmcl.compare_homolog_groups --cov-report html -n 8 --cov-config .coveragerc -p no:cacheprovider -p no:logging --durations=10 -m $@
