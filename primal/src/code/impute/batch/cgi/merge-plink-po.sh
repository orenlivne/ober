#!/bin/bash
#-----------------------------------------------------------------
# Merge all parent-of-origin imputed plink files into a master
# bed file.
#-----------------------------------------------------------------

# Input arguments
imputed_data_dir="$1" # Directory of final imputed data and PO data

OUT_DIR="${OBER_OUT}/impute_cgi/data-sets"
PLINK="plink --noweb --nonfounders"

#-----------------------------------------------------------------
# Merge all chromosomes into one big PLINK data set (one for
# regular genotypes, one for PO genotypes)
#-----------------------------------------------------------------
if false; then
cd ${imputed_data_dir}
for t in imputed_cgi imputed_cgi.po; do
    rm -f merge-list
    for c in `seq 2 22`; do
	echo "$t.chr$c.bed $t.chr$c.bim $t.chr$c.fam" >> merge-list
    done
    ${PLINK} --bfile $t.chr1 --merge-list merge-list --make-bed --out $t
    rm -f merge-list
done
fi

#-----------------------------------------------------------------
# Create subset of LD-pruned SNPs Catherine's PO GWAS
#-----------------------------------------------------------------
mkdir -p ${OUT_DIR}/qc-pruned-po
cd ${OUT_DIR}/qc-pruned-po
awk '{print $2}' ${OUT_DIR}/qc-pruned/qc.pruned.bim > qc-pruned.snplist

${PLINK} --bfile ${imputed_data_dir}/imputed_cgi.po --extract qc-pruned.snplist --make-bed --out qc-pruned-po
