#!/bin/bash
#-----------------------------------------------------------------
# Run impute2 genotype QC on a single node.
# Usage: qc-impute2.sh <input-chunk> <min-hets> <num-bins>
#-----------------------------------------------------------------

cat $1 | parallel -j 24 qc-impute2 -n $2 -b $3 -p
