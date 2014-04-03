#!/usr/bin/env python
'''
============================================================
Read a PLINK BIM file and compute part size per chromosome.

Split a PLINK BED file into chromosome equal bp-length parts
accordingly.

Recode alleles so that 1=minor, 2=major.

Prerequisites: PLINK must be installed and be on the PATH.

Created on Septmeber 22, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, csv, os, shutil, traceback, math, numpy as np, util, impute.batch.batch_util as bu, \
impute.data.constants as constants
from optparse import OptionParser
from impute.data.constants import MEGA_BASE_PAIR
from util import mkdir_if_not_exists

#---------------------------------------------
# Constants
#---------------------------------------------
# Program name
PROGRAM = os.path.basename(sys.argv[0])

#---------------------------------------------
# Methods
#---------------------------------------------

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __chr_endpoints(bim_file):
    '''Return a dictionary of chromosome-to-bp-length from a BIM file. Note: the relevant length
    is measured between the first and last SNPs, and may be shorter than the total chromosome length.'''
    endpoints = {}
    for line in csv.reader(bim_file, delimiter='\t', skipinitialspace=True):
        if line:
            (key, bp) = (int(line[0]), int(line[3]))
            value = endpoints.setdefault(key, [np.inf, -1])
            if bp < value[0]:
                value[0] = bp
            if bp > value[1]:
                value[1] = bp
    # Include only relevant chromosomes (CHROMOSOMES)
    return dict((k, endpoints[k]) for k in constants.CHROMOSOMES)

def __chr_length(endpoints):
    '''Convert chromosome endpoint dictionary to length. Unit is bp.'''
    return dict((k, v[1] - v[0]) for (k, v) in endpoints.iteritems())

def __chr_part_count(key_count, part_size):
    '''Convert key counts to part counts.'''
    return dict((k, __part_count(v, part_size)) for (k, v) in key_count.iteritems())

def __part_count(num_lines, part_size):
    return int(math.ceil((1.0 * num_lines / part_size)))

def __write_int_dict(part_count, out):
    '''Write part counts per chromosome to out_file.'''
    for (key, part) in part_count.iteritems():
        out.write('%d %d\n' % (key, part))
        
def __split_chromosome(in_file, chrom, part_size, fout):
    '''Split a single chromosome's data into num_part files of at most part_size lines each.'''
    # print 'Splitting chromosome to out file', [out.name for out in fout]
    (part_num, count) = (0, 0)
    for line in csv.reader(open(in_file, 'rb'), delimiter=' ', skipinitialspace=True):
        if line:
            (key, value) = (int(line[0]), int(line[1]))
            if key == chrom:
                fout[part_num].write('%d %d\n' % (key, value))
                count += 1
                if count == part_size:
                    count = 0
                    part_num += 1 

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    usage = 'Usage: %s <plink_set> <part_size>\n\n' \
        'Read plink_set.bim file and compute part size per chromosome. Save in plink_set.part;\n' \
        'or split a PLINK BED file into chromosome parts accordingly.\n' \
        'part_size is the slice size in Mb.\n\n' \
        'Example 1: Print part counts by chromosome: %s -s genome 2\n' \
        'Example 2: Split chromosome #2 into 2-Mb parts: %s -c 2 genome 2\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM, PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-v', '--debug'        , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-c', '--chr'          , type='int'           , dest='chrom', default=0,
                      help='Chromosome number to split')
    parser.add_option('-s', '--summary'      , action='store_true'  , dest='summary', default=False,
                      help='Print part counts by chromosome')
    parser.add_option('-r', '--recode'      , action='store_true'  , dest='recode', default=True,
                      help='Recode alleles to 1=minor, 2=major (if False, a random assignment to 1,2 is made)')
    parser.add_option('-o', '--out'          , type='str'           , dest='out_base_name', default=None,
                      help='Output PLINK data set base name')
    parser.add_option('-g', '--out-gxn'     , type='str'           , dest='out_gxn', default=bu.ARG_NONE,
                      help='Output directory of GXN files (if specified, copies FRQ files to that directory)')
    (options, args) = parser.parse_args(sys.argv[1:])
    if options.out_gxn.startswith(bu.ARG_NONE):
        options.out_gxn = ''
    #print 'Using target dir out_gxn = ''%s''' % (options.out_gxn,) 
    if not (options.summary ^ (options.chrom != 0)):
        print usage
        print('\nMust specify either summary mode or chromosome number')
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
    if len(args) != 2:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    try:
        (base_name, part_size) = (args[0], int(args[1]) * MEGA_BASE_PAIR)
        if not options.out_base_name:
            options.out_base_name = base_name
        mkdir_if_not_exists(os.path.dirname(options.out_base_name))
        
        endpoints = __chr_endpoints(open(base_name + '.bim', 'rb'))
        part_count = __chr_part_count(__chr_length(endpoints), part_size)
        if options.summary:
            # Compute chromosome length in BP, compute part counts, save them to file/stdout        
            with (open(bu.partcountname(options.out_base_name), 'wb')) as out:
                __write_int_dict(part_count, out)
            np.savez(bu.endpointcountname(options.out_base_name),
                     endpoints=np.array([(x, y) for (k, (x, y)) in endpoints.iteritems()]), 
                     part_size=np.array([part_size]))
        else:
            # Extract specific chromosome into output file
            part_names = bu.partnames(bu.chrnames(options.out_base_name, part=options.chrom),
                                      num_parts=part_count[options.chrom])
            part_names_gxn = bu.partnames(bu.chrnames(options.out_gxn, part=options.chrom),
                                          num_parts=part_count[options.chrom])
            (start, stop) = endpoints[options.chrom]
            endpoints = util.brange(start, stop, part_size, endpoint=True, dtype=int)
            for (part, part_start) in enumerate(endpoints[:-1]):
                out = part_names[part]
                plink_cmd_base = '%s --cm --bfile %s --chr %d --from-bp %s --to-bp %s --out %s' \
                % (bu.PLINK, base_name, options.chrom, part_start, endpoints[part + 1], out)
                
                if options.recode:
                    # First, compute allele frequencies with PLINK  
                    util.run_command(plink_cmd_base + ' --nonfounders --freq')
                    # Convert frequencies file that to a reference allele recoding
                    # file (a file containing the list of SNPs and their minor allele letter)
                    bu.frq_to_minor_file(out + '.frq', out + '.mnr') 
                    if options.out_gxn:
                        # Copy FRQ to target output directory
                        out_frq = part_names_gxn[part] + '.frq'
                        mkdir_if_not_exists(os.path.dirname(out_frq))
                        shutil.copy(out + '.frq', out_frq)
    
                    # Then convert binary PLINK to a 12-recoded TPED, where 1=minor allele for each SNP                 
                    cmd = '%s --transpose --recode12 --reference-allele %s.mnr' % (plink_cmd_base, out)
                    util.run_command(cmd)
                else:
                    # No recoding, just convert binary to 2-recoded TPED. PLINK assigns "1" to
                    # the first-encountered allele in the file for each SNP.
                    cmd = '%s --transpose --recode12' % (plink_cmd_base,)
                    util.run_command(cmd)                    
                
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
