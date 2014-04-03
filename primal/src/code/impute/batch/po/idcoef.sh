#!/bin/bash
#-------------------------------------------------------------------------
# Calculate detailed identity coefficients from a haplotype IBD
# segment index.
#
# Author: Oren E. Livne
# Date:   01-DEC-2012
#-------------------------------------------------------------------------

#-----------------
# Input parameters
#-----------------
# Local imputed data output directory
segment_index="$1" # Location of segment index
out="$2" # Output directory
start_chr="$3"
stop_chr="$4"
num_jobs="$5" # #jobs to run in parallel

# Calculate identity state counts for a single chromosome
function do_chrom
{
    # Location of segment index
    local segment_index="$1"
    # Output directory
    local out="$2"
    # Chromosome number
    local chrom="$3"

    echo "Chromosome ${chrom}"
    ${OBER_CODE}/impute/impute/poo/idcoef.py -c ${chrom} ${segment_index} ${out}/idcoefs-chr${chrom}.txt
}

#---------------------
# Main program
#---------------------
mkdir -p ${out}
export -f do_chrom
monitor-usage 30 false >& ${out}/monitor-usage.out &
if [ ${num_jobs} -eq 1 ]; then
    # Serial run
    for chrom in `seq ${start_chr} ${stop_chr}`; do
	do_chrom ${segment_index} ${out} ${chrom}
    done
else
    # Parallel run
    seq ${start_chr} ${stop_chr} | parallel -j ${num_jobs} do_chrom ${segment_index} ${out}
fi
