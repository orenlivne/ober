#!/bin/bash
#-----------------------------------------------------------------
# Remove trailing tab in CGI output files.
#
# Replace empty fields with '-'. Assumes there are no more than
# 2 consecutive empty fields.
# 
# Removes partial calls per sequencing meeting decision to use
# only genotypes called as HH by CGI.
# William: "Alleles called as VQLOW by CGI are indicated as
# no-calls in Jessica's file [all....testvar...tsv.gz]. E.g.:
#
# HH -> 00, 01, 11
# HL -> 0N, 1N
# LL -> NN, NN."
#-----------------------------------------------------------------
if [ $# -ge 1 ]; then
    start_chrom=$1
    stop_chrom=$1
else
    start_chrom=1
    stop_chrom=22
fi

out=/lustre/beagle/ober/users/oren/out/impute_cgi/genotypes
mkdir -p ${out}
for c in `seq ${start_chrom} ${stop_chrom}`; do
#   sed 's/.\{1\}$//' all/imputed_cgi.chr$c.tsv > trimmed/imputed_cgi.chr$c.tsv &
#    sed -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' all/imputed_cgi.chr$c.tsv > trimmed/imputed_cgi.chr$c.tsv &
    zcat $OBER_DATA/cgi/all.2012-09-20.testvar.chr$c.tsv.gz | sed -e 's/\t\t/\t-\t/g' -e 's/\t\t/\t-\t/g' -e 's/\t0N/\tNN/g' -e 's/\t1N/\tNN/g' > ${out}/genotypes.chr$c.tsv &
done
wait
