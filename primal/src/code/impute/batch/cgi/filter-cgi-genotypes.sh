#!/bin/bash
#-----------------------------------------------------------------
# Filter CGI genotypes to have at least threshold call rate.
#-----------------------------------------------------------------
chrom="$1"
threshold="$2"

exec="${OBER_CODE}/impute/batch/cgi/filter-cgi-genotypes"
in_prefix="$OBER_DATA/cgi/all.2012-09-20.testvar"
out_dir="$OBER_OUT/impute_cgi/filtered-genotypes/${threshold}"

${exec} ${in_prefix}.chr${chrom}.tsv.gz ${threshold} > ${out_dir}/variant_id.chr${chrom}
