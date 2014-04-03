#!/bin/bash
#--------------------------------------------------------------------
# Count genotypes in an imputed CGI file of a single chromosome.
#
# Author: Oren E. Livne
# Date:   01-DEC-2012
#--------------------------------------------------------------------
chrom="$1"
out="$2"
opts="$3"

cd ${OBER_OUT}/impute_cgi/final
out_file=${out}.chr${chrom}.txt
touch ${out_file}
cat imputed_cgi.chr${chrom}.tsv | python ${OBER_CODE}/impute/impute/cgi/count_cgi_genotypes.py ${opts} > ${out_file}
