#!/bin/bash
#-----------------------------------------------------------------
# Split IBD segments, to be fed into the segment indexing process.
#
# Usage: split-segments.sh <chrom> <nodes> <instances_per_node> <region_size> <#splitting jobs to run in parallel>
#-----------------------------------------------------------------
chrom="$1"
nodes="$2"
instances_per_node="$3"
region_size="$4"
split_jobs="$5"
work="$6"

# Copy input file to RAMDISK
ramdisk="/dev/shm"
rm -f $ramdisk/segments.out >& /dev/null
cp $OBER_OUT/phasing/chr$chrom/segments.out $ramdisk

# Run split segments
echo $OBER_CODE/impute/batch/index_segments/split-segments -j $split_jobs -a 4 -b 4 $chrom $ramdisk/segments.out $work/chr$chrom/index_segments/node- "" "index_segments-" $nodes $instances_per_node $region_size
$OBER_CODE/impute/batch/index_segments/split-segments -j $split_jobs -a 4 -b 4 $chrom $ramdisk/segments.out $work/chr$chrom/index_segments/node- "" "index_segments-" $nodes $instances_per_node $region_size

# Clean RAMDISK
rm -f $ramdisk/segments.out >& /dev/null
