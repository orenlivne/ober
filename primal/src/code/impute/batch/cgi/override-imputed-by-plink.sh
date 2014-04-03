#!/bin/bash
#--------------------------------------------------------------------
# Override imputed genotypes by genotypes read from a PLINK BED file.
# For instance, this could be the Affymetrix SNPs.
#
# If a TPED file exists (assumed to contain the chromosome of
# interest's data only), use it directly instead of the BED file.
#
# Author: Oren E. Livne
# Date:   09-MAY-2013
#--------------------------------------------------------------------
# Constants
EXEC="${OBER_CODE}/impute/impute/cgi/override_imputed_by_plink.py -b 20000"
out_temp="/dev/shm/imputed-override"

# Read input arguments
DARGS=65
PROGNAME=`basename $0`
if [[ ( $# -lt 5 ) || ( $# -gt 6 ) ]]; then
  echo "Usage: ${PROGNAME} <in-file-prefix> <out-dir> <plink_prefix> <chrom> <temp-dir> [override_flags]"
  echo ""
  echo "Override imputed genotypes by genotypes read from a PLINK file. If chromosome is empty,"
  echo "processes all chromosomes in parallel."
  exit $E_BADARGS
fi
in_prefix="$1"
out="$2"
plink_prefix="$3"
chrom="$4"
out_temp="$5"
flags=""
if [[ $# -ge 6 ]]; then
    flags="$6"
fi

rm -rf ${out_temp}
mkdir -p ${out} ${out_temp}

# Override a single chromosome file.
function override_chrom
{
    local chrom="$1"
    suffix="chr${chrom}.tsv"
    cd ${out_temp}

    if [ ! -f imputed_cgi.${suffix}.gz.tbi ]; then
        # Copy data into temp directory (under ramdisk)
	cp ${in_prefix}.${suffix}.gz ${out_temp}

        # Unzip the imputed data file
	bgzip -d imputed_cgi.${suffix}.gz
	
	if [ ! -f ${plink_prefix}.tped ]; then
            # Extract this chromosome's SNP from the plink file
	    plink --bfile ${plink_prefix} --noweb --recode12 --transpose --chr ${chrom} --out plink.${chrom}
	    
            # Create overridden-imputed file
            ${EXEC} ${flags} plink.${chrom} imputed_cgi.${suffix} > imputed_override.${suffix}
	else
	    # TPED exists ==> assuming TPED is already made for this chromosome; use it
	    ${EXEC} ${flags} ${plink_prefix} imputed_cgi.${suffix} > imputed_override.${suffix}
	fi

        # Index overridden-imputed file
	bgzip -c imputed_override.${suffix} > imputed_cgi.${suffix}.gz
	tabix -s 2 -b 3 -e 4 imputed_cgi.${suffix}.gz
    fi

    # Move output to final destination
    #rm -f ${out}/imputed_cgi.${suffix}.gz*
    cp imputed_cgi.${suffix}.gz* ${out}

    # Clean temporary files
#    rm -rf ${out_temp}
}

# Main program

if [[ "x${chrom}" != "x" ]]; then
    override_chrom ${chrom}
else
    for chrom in `seq 1 22`; do
	override_chrom ${chrom} &
    done
    wait
fi