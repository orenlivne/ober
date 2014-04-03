#!/bin/bash
#-------------------------------------------------------------------------
# Extract a +-1Mb window of imputed PO data around each of several GWAS
# hits.
#
# Input parameters:                                                                                                   
# $1 = in_file = list of hits (CGI variant IDs)
# $2 = out_dir = output directory
# $3 = num_jobs = number of jobs to run in parallel 
#
# Author: Oren E. Livne
# Date:   13-MAR-2014
#-------------------------------------------------------------------------
# note:  this program can be be made much faster by extracting
# from chromosomal files rather than from the master PLINK file.
#-------------------------------------------------------------------------

#-----------------
# Input parameters
#-----------------
in_file="$1"
out_dir="$2"
num_jobs="$3"

# Calculate identity state counts for a single chromosome
function do_hit
{
    master_imputed_file="${OBER_OUT}/impute_cgi/imputed-override3/imputed_cgi.po"
    out_dir="$1"
    hit="$2"
    echo "Running plink --noweb --nonfounders --bfile ${master_imputed_file} --snp ${hit} --window 1000 --out ${out_dir}/imputed_cgi.po.${hit} --make-bed"
    plink --noweb --nonfounders --bfile ${master_imputed_file} --snp ${hit} --window 1000 --out ${out_dir}/imputed_cgi.po.${hit} --make-bed
}

#---------------------
# Main program
#---------------------
mkdir -p ${out_dir}
export -f do_hit
if [ ${num_jobs} -eq 1 ]; then
    # Serial run
    cat ${in_file} | ( while read hit; do do_hit ${out_dir} ${hit}; done ) 
else
    # Parallel run
    cat ${in_file} | parallel -j ${num_jobs} do_hit ${out_dir}
fi
