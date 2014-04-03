'''
============================================================
Centralizes global constants.

Created on May 31, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, csv, numpy as np, util
from collections import OrderedDict
from impute.data.constants import CHROMOSOMES

#---------------------------------------------
# Constants
#---------------------------------------------
'''Types of data set partitions.'''
PART, CHROMOSOME = range(2)
# File name part separator
SEPARATOR = '_'

'''PLINK commands'''
# Base command
PLINK = 'plink --noweb' 

# A command-line argument prefix that indicates no value
ARG_NONE = 'NONE'

#---------------------------------------------
# File Names
#---------------------------------------------
def partcountname(base_name):
    '''Given a base file name, generate the standard part count file name.'''
    (prefix, _) = os.path.splitext(base_name)
    return prefix + '.part'

def endpointcountname(base_name):
    '''Given a base file name, generate the standard chromosome bp endpoint file name.'''
    (prefix, _) = os.path.splitext(base_name)
    return prefix + '.range.npz'

def chrnames(base_name, parts=CHROMOSOMES, part=None):
    '''Given a base file name, generate the standard part file names.  If chromosome is not None,
    outputs a single chromosome name.'''
    return __partnames(base_name, parts, CHROMOSOME, part=part)

def partnames(base_name, num_parts=None, parts=None, part=None):
    '''Given a base file name, generate the standard chromosome part file names. If part is not None,
    outputs a single part name.'''
    return __partnames(base_name, parts if num_parts is None else xrange(0, num_parts), PART, part=part)

def partition_names(base_name, part_type, parts=None, part=None):
    '''Generalizes chrnames, partnames.'''
    if part_type == CHROMOSOME:
        return chrnames(base_name, parts=parts, part=part)
    elif part_type == PART:
        return partnames(base_name, parts=parts, part=part)
    else:
        raise ValueError('Unsupported partition type ''%s''' % (part_type,))
    
def loadparts(file_name):
    '''Load a chromosome-to-part dictionary from a text file.'''
    part_num = np.loadtxt(file_name).astype('int')
    return dict(zip(part_num[:, 0], part_num[:, 1]))

#---------------------------------------------
# Other Public Methods
#---------------------------------------------
def frq_to_minor_file(frq_file_name, minor_file_name):
    '''Convert a PLINK frequency file to a PLINK minor allele coding file.'''
    with open(minor_file_name, 'wb') as out:
        for line in csv.reader(open(frq_file_name, 'rb'), delimiter=' ', skipinitialspace=True):
            if line: out.write('%s %s\n' % (line[1], line[2]))

def get_bim_metadata(bim_file):
    '''Return (a)
    a dictionary of chromosome-to-bp-length from a BIM file. Note: the relevant length
    is measured between the first and last SNPs, and may be shorter than the total chromosome length;
    (b) a dictionary of chromosome-to-list-of-snp-names.'''
#    a = list(itemutil.index_of_change(((int(x[0]), int(x[3])) 
#                                       for x in csv.reader(bim_file, delimiter='\t', skipinitialspace=True)), 
#                                      output_first=True, output_last=True, 
#                                      output_value=True, 
#                                      comparator=lambda x,y: x[0]==y[0]))
#    endpoints = np.diff(np.array([(x[1][1] if x[1] else 0, x[2][1] if x[2] else 0) for x in a]).flatten())[1:2*22:2]

    # A BIM file is sorted by chromosome, then by SNP base-pair location
    snp_names = util.mdict()
    chr_endpoints = {}
    prev = start = None
    for line in csv.reader(bim_file, delimiter='\t', skipinitialspace=True):
        chrom = int(line[0])
        curr = (chrom, int(line[3]))
        snp_names[chrom] = line[1]
        if not start: start = curr
        if prev and chrom != prev[0]:
            chr_endpoints[start[0]] = prev[1] - start[1]
            start = curr
        prev = curr
    chr_endpoints[start[0]] = prev[1] - start[1] # Last chromosome
    return dict((k, chr_endpoints[k]) for k in CHROMOSOMES if chr_endpoints.has_key(k)), snp_names

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __partnames(base_name, parts, part_type, part=None):
    '''Given a base file name, generate the standard part file names. If part is not None,
    outputs a single part name.'''
    (prefix, suffix) = os.path.splitext(base_name)
    suffix_fmt = __suffix_fmt(part_type)
    name = lambda x: (prefix + suffix_fmt + suffix) % (x,) 
    return OrderedDict((x, name(x)) for x in parts) if part is None else name(part)

def __suffix_fmt(part_type):
    if part_type == CHROMOSOME:
        return SEPARATOR + 'chr%d'
    elif part_type == PART:
        return SEPARATOR + 'part%d'
    else:
        raise ValueError('Unsupported partition type ''%s''' % (part_type,))
