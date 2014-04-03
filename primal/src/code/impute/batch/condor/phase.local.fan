#-----------------------------------------------------------
# Phasing condor_fan pipeline
# Suitable for a local/shared file system across the cluster
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
# Fill missing genotypes (0=none, 1=impute, 2=impute & fill rest with random)
impute=%(arg7)s
# Phasing stage (0 = real run; 1-5 = run only this phasing stage)
stage=%(arg8)s
# Gaixin-format output directory
out_gxn=%(arg9)s
batch_home=%(arg10)s

# Optional global place-holders - Convenient aliases
input_path=%(input_dir)s/%(input)s
[END DEFAULT] 					# Do not delete this line

[preprocess]
executable = none

[split_chr]
executable = %(python)s
arguments = %(batch_home)s/split.py -g %(out_gxn)s/chr$(chr)/%(input)s_phased_chr$(chr) -c $(chr) %(input_path)s %(part_size)s -o chr$(chr)/%(input)s

[map]
executable = %(python)s
arguments = %(batch_home)s/map.py -g %(out_gxn)s/chr$(chr)/%(input)s_phased_chr$(chr)_part$(part) -f %(impute)s ../split_chr/chr$(chr)/%(input)s_chr$(chr)_part$(part) %(input_path)s.pdg.tfam chr$(chr)/%(input)s_phased_chr$(chr)_part$(part) -v -t %(stage)s

[reduce_chr]
executable = %(python)s
arguments = %(batch_home)s/reduce.py -s 0 -e $(num_parts) ../map/chr$(chr)/%(input)s_phased_chr$(chr) part chr/%(input)s_phased_chr$(chr)

[reduce]
executable = %(python)s
arguments = %(batch_home)s/reduce.py -s %(start_chr)s -e %(stop_chr_p1)s ../reduce_chr/chr/%(input)s_phased chr result/%(input)s_phased
