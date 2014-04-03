#!/bin/bash
#-----------------------------------------------------------------
# Create common non-RS variants data set for Catherine's GWAS.
#
# Assumes the list of variants is in
# ${OUT_DIR}/common-novel-qc.txt
#-----------------------------------------------------------------

OUT_DIR="$1"
PLINK="plink --noweb --nonfounders"

mkdir -p ${OUT_DIR}
cd ${OUT_DIR}
awk '{print $1}' common-novel-qc.txt > common-novel-qc.snplist
${PLINK} --bfile ${im2}/imputed_cgi --extract common-novel-qc.snplist --make-bed --out common-novel-qc
