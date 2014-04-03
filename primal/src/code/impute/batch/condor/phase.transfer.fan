#-----------------------------------------------------------
# Phasing condor_fan pipeline
# Uses file transfer
#-----------------------------------------------------------
# Pre-defined variables:
# pwd - current working directory 

[DEFAULT]
# Mandatory place-holders
# PLINK input data set base name
# Number of times to retry a processing job before failing
num_retrys=2
# A file containing chromosome-to-#parts mapping
part_count=%(arg4)s

# Optional global place-holders - phasing parameters
# in Mb
part_size=50
# Phasing stage (0 = real run; 1-5 = run only this phasing stage)
stage=1
# Python executable
python = %(arg0)s
# Input data set
input_dir=%(arg1)s
input = %(arg2)s
# Output directory
out=%(arg3)s
# Start chromosome
start_chr=%(arg5)s
stop_chr_p1=%(arg6)s
# Fill missing genotypes (map.py flag; string should be "-f" if on, "" if off)
fill_missing=%(arg7)s
# Phasing stage (0 = real run; 1-5 = run only this phasing stage)
stage=%(arg8)s

# Optional global place-holders - Convenient aliases
input_path=%(input_dir)s/%(input)s
[END DEFAULT] 					# Do not delete this line

[preprocess]
executable = none

[split_chr]
executable = %(python)s
arguments = split.py -c $(chr) %(input_path)s %(part_size)s -o chr$(chr)/%(input)s
transfer_input_files = %(pwd)s/split.py, %(input_path)s.bed, %(input_path)s.bim, %(input_path)s.fam, %(input_path)s.pdg.tfam
transfer_output_files = chr$(chr)

[map]
executable = %(python)s
arguments = map.py %(fill_missing)s %(input)s_chr$(chr)_part$(part) %(input)s.pdg.tfam chr$(chr)/%(input)s_phased_chr$(chr)_part$(part) -v -s %(stage)s
transfer_input_files = %(pwd)s/map.py, ../split_chr/chr$(chr)/%(input)s_chr$(chr)_part$(part).tped, ../split_chr/chr$(chr)/%(input)s_chr$(chr)_part$(part).tfam, %(input_path)s.pdg.tfam 
transfer_output_files = chr$(chr)

[reduce_chr]
executable = %(python)s
arguments = reduce.py -s 0 -e $(num_parts) chr$(chr)/%(input)s_phased_chr$(chr) part chr/%(input)s_phased_chr$(chr)
transfer_input_files = %(pwd)s/reduce.py, ../map/chr$(chr)
transfer_output_files = chr

[reduce]
executable = %(python)s
arguments = reduce.py -s %(start_chr)s -e %(stop_chr_p1)s chr/%(input)s_phased chr result/%(input)s_phased
transfer_input_files = %(pwd)s/reduce.py, ../reduce_chr/chr
transfer_output_files = result
