#!/bin/bash
#-----------------------------------------------------------------
# Merge all imputed plink files into a master bed file.
#-----------------------------------------------------------------

OUT_DIR="${OBER_OUT}/impute_cgi/data-sets"
PLINK="plink --noweb --nonfounders"

#-----------------------------------------------------------------
# Merge all chromosomes into one big imputed data set under ${im2}
#-----------------------------------------------------------------
cd ${im2}
rm -f merge-list
for c in `seq 2 22`; do
    echo "imputed_cgi.chr$c.bed imputed_cgi.chr$c.bim imputed_cgi.chr$c.fam" >> merge-list
done
${PLINK} --bfile imputed_cgi.chr1 --merge-list merge-list --make-bed --out imputed_cgi

#-----------------------------------------------------------------
# Create common non-RS variants data set for Catherine's GWAS
#-----------------------------------------------------------------
${OBER_CODE}/impute/batch/cgi/common-novel-qc.sh ${OUT_DIR}/common-novel-qc

#-----------------------------------------------------------------
# Create data set with all variants that passed QC for Gorka
#-----------------------------------------------------------------
mkdir -p ${OUT_DIR}/qc
cd ${OUT_DIR}/qc
${PLINK} --bfile ${im2}/imputed_cgi --extract qc.txt --make-bed --out qc
# Split to chromosomal files per his request
for c in `seq 1 22`; do
    ${PLINK} --bfile qc --out qc.chr$c --make-bed --chr $c
done
# Calculate frequencies
${PLINK} --bfile qc --out qc --freq
