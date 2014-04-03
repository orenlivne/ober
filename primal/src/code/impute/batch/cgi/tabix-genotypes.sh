#!/bin/bash
#-----------------------------------------------------------------
# Compress and index genotype files using tabix.
# Input arguments:
# $1: file name pattern. Format: input_dir/file.chr#chrom#.tsv
#     the strong #chrom# is successively replaced by 1..22.
#-----------------------------------------------------------------
name="$1"
for chrom in `seq 1 22`; do
  file=`echo $name | sed "s/#chrom#/${chrom}/g"`;
  echo "Tabixing $file ..."
  ( bgzip -c $file > $file.gz ; tabix -s 2 -b 3 -e 4 $file.gz ) &
done
wait
