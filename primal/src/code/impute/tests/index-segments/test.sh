#!/bin/bash
# Test segment indexing at a problematic monogenic variant (closest SNP: chr18:index 2609)

TEST_DIR="${OBER_CODE}/impute/tests/index-segments"
EXEC="${OBER_CODE}/impute/impute/ibd/index/index_segments_beagle.py"
PHASED="${OBER_OUT}/phasing/chr18/hutt.phased.info.npz"
SEGMENTS="${TEST_DIR}/index-segments-chr18-2600-2700.out" #"${TEST_DIR}/test-0000-100-samples.dat"

python ${EXEC} -s 2609 -w 0 -f -p 1 -v 2 -a amg -l 0.4 -r 100 ${TEST_DIR}/test-2.in ${PHASED} ${SEGMENTS} ${TEST_DIR}
#python ${EXEC} -f -p 1 -v 1 -a amg -l 0.4 -r 100 ${TEST_DIR}/test-2.in ${PHASED} ${SEGMENTS} ${TEST_DIR}
