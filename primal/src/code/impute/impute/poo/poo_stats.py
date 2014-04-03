#!/usr/bin/env python
'''
============================================================
Generate some POO statistics for a single variant's
imputed genotype data file. 

Created on January 22, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, os, sys, itertools as it, matplotlib.pyplot as P
from optparse import OptionParser
from collections import Counter
from statutil import cumbin
from impute.data.io_pedigree import read_sample_id

#---------------------------------------------
# Constants
#---------------------------------------------
'''# of metadata columns in an itabix line.'''
NUM_METADATA_COLS = 8 

#---------------------------------------------
# Methods
#---------------------------------------------
TOTAL, QF, NON_QF = xrange(3)

####################################################################################
class PooStats(object):
    '''Generate POO statistics per each SNP in an itabix output file.'''
    __CGI_LETTER_TO_ALLELE = {'N': '0', '0': '1', '1': '2'}  # String-to-string coding, slightly more convenient here

    def __init__(self, file_name, affy_bim, genotyped_id_file, debug=False): 
        p = im.hutt_pedigree()
        self.qf = p.quasi_founders
        self.non_qf = np.setdiff1d(xrange(p.num_genotyped), p.quasi_founders)
        self.debug = debug
        # Load affy SNP names 
        self.affy = set(np.loadtxt(affy_bim, usecols=[1], dtype=str))
        self.sample_id = read_sample_id(genotyped_id_file)
        # Load data, cache statistics in the data field
        self.data = self.__stats_struct(file_name)
    
    def print_stats(self, out=sys.stdout):
        '''Print a SNP statistics report to stdin.'''
        # Print header line
        out.write('%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
        ('#', 'Chrom', 'Position', 'Name', 'Affy?', 'CR', 'MAF', 'Total Counts', 'QF Counts', 'Non-QF Counts'))
        out.write('%s\t%s\t\t%s\t%s\t\t\t%s\t%s\t%s\t%s\t%s\n' % \
        ('', '', '', '', '', '', '10:01 (no-PO)', '10:01 (no-PO)', '10:01 (no-PO)'))
        out.write('-' * 125 + '\n')
        for i, s in enumerate(self.data):
            chrom, bp, name, is_affy, call_rate, maf, count = s        
            out.write('%-2d\t%-2d\t%-9d\t%-10s\t%-2s\t%-5.3f\t%-5.3f\t%3d:%3d (%3d) %.2e\t%3d:%3d (%3d)\t%3d:%3d (%3d)\n' % \
            (i, chrom, bp, name, yesno(is_affy), call_rate, maf,
             count[0][0], count[0][1], count[0][2], p_value(count[0][0], count[0][1]),
             count[1][0], count[1][1], count[1][2],
             count[2][0], count[2][1], count[2][2]))

    def __stats(self, file_name):
        '''A generator of POO statistics per SNP loaded from the tabix file file_name. Streams the data.'''
        allele = PooStats.__CGI_LETTER_TO_ALLELE
        for i, line in enumerate(x.rstrip().split() for x in open(file_name, 'rb')):
            if self.debug and i % 1000 == 0: print 'Loaded %d SNPs' % (i,)
            # Read metadata
            chrom, name, bp = int(line[1][3:]), PooStats.__parse_rs_number(line[7]), int(line[3])
            # Read imputed genotypes. Each is a three-character string: tag,a1,a2
            g = np.array([x[0] + allele[x[1]] + allele[x[2]] for x in line[8:]])
            
            # Calculate genotype statistics
            is_affy = name in self.affy
            call_rate = float(len(np.where(['N' not in x for x in g])[0])) / len(g)
            count = [(c['221'], c['212'], c['012'] + c['021'] + c['112'] + c['121']) for c in
                      (Counter(genotypes) for genotypes in (g, g[self.qf], g[self.non_qf]))]
            
            # Allele frequencies
            a = Counter(x[i] for x in g for i in xrange(1, 3))
            a1, a2 = a['1'], a['2']
            maf = -1 if a1 == 0 and a2 == 0 else float(min(a1, a2)) / (a1 + a2)
            
            yield chrom, bp, name, is_affy, call_rate, maf, count

    def __stats_struct(self, file_name):
        '''Loads all stats() into a record array.'''
        return np.array(list(self.__stats(file_name)),
                        dtype=[ ('chrom', 'i2'),  # Chromosome number
                                ('bp', 'i10'),  # Base-pair position
                                ('name', 'S12'),  # SNP name
                                ('is_affy', np.bool),  # Is SNP on affy array?
                                ('call_rate', np.float),  # Imputation call rate
                                ('maf', np.float),  # Minor allale frequency
                                ('count', '(3,3)i4'),  # Het counts (total, Quasi-Founder (QF) hets, non- QF hets; for each, 201, 210, missing-PO-hets)
                                ])

    def p_value(self, count_type):
        '''Return the two-tailed binomial p-value of the 10:01 ratio for
        count_type=0 (total), 1 (QF) or 2 (non-QF).'''
        c = self.data['count']
        return np.array([p_value(x, y) for x, y in it.izip(c[:, count_type, 0], c[:, count_type, 1])])

    def scatter_pvalue(self, count_type, min_hets=10):
        '''Plot -log10(p-value) along the chromosome. count_type=0 (total), 1 (QF) or 2 (non-QF).
        Include variants with at least min_hets hets.'''
        c = self.data['count']
        called, missing = c[:, count_type, 0] + c[:, count_type, 1], c[:, count_type, 2]
        i = np.where(called >= min_hets)
        called, missing, p, bp = called[i], missing[i], self.p_value(count_type)[i], self.data['bp'][i]
        # Scale luminance (<=> index into colormap) to [0,1]
        # luminance = np.digitize(called.astype(float) / (called + missing), [0, 0.2, 0.5, 0.8, 0.95, 1])
        # luminance = 1 - luminance / max(luminance)
        luminance = called.astype(float) / (called + missing)
        # luminance = luminance ** 4  # Space out higher luminance values 
        cm = P.cm.get_cmap('Spectral')
        h = P.scatter(bp, -np.log10(p), lw=0, s=30, color='b', cmap=cm, c=luminance)
        P.colorbar(h)
        P.xlabel('Base Pair')
        P.ylabel('log10(p-value)')
        P.title('10:01 2-Tailed Significance along Chromosome %d' % (self.data['chrom'][0],))
        P.ylim(-0.5, P.ylim()[1])
        P.xlim(min(bp) - 1e6, max(bp) + 1e6)

    @staticmethod
    def __parse_rs_number(s):
        '''Parse a CGI xref column value into the correspodning latest dbSNP RS number.''' 
        return '-' if s == '-' else max(((int(x[0][6:]), x[1]) if x[0].startswith('dbsnp') else (-1, None)) for x in (x.split(':') for x in s.split(';')))[1]

'''Convert a boolean to a Y/N string.'''
yesno = lambda x: 'Y' if x else 'N'

def p_value(c1, c2):
    '''Significance of observing c1:c2 ratio. One-tailed.'''
    if c1 + c2 == 0: return 1
    return 1 - cumbin(c1 + c2, 0.5, max(c1, c2)) + cumbin(c1 + c2, 0.5, min(c1, c2))

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <genotype-file>\n' \
        'Generate parent-of-origin statistics for a single imputed variant.\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
#    parser.add_option('-m', '--flag-maf-not-one', action='store_true'  , dest='flag_maf_not_one', default=False,
#                      help='If true, list only the SNPs in which the 1-allele is not the minor one')
    options, args = parser.parse_args(sys.argv[1:])
    genotype_file_name = args[0]
    if len(args) != 1:
        print usage
        sys.exit(1)
    
    # Load genotype data
    print 'File: %s' % (args[0],)
    s = PooStats(genotype_file_name, os.environ['OBER_DATA'] + '/hutt/hutt.3chipoverlap.clean.bim', os.environ['OBER_DATA'] + '/hutt/hutt.3chipoverlap.clean.fam')

    P.figure(1)
    P.clf()
    s.scatter_pvalue(TOTAL)
    P.show()
    # P.savefig('/home/oren/ober/po-p.png')
