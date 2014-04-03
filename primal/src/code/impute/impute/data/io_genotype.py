'''
============================================================
Read and write genotype sets from/to a file set. Currently
supports only the PLINK format.

Requires:
- tped file
- tfam file

File format - see http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
     
Example:
TFAM file:
HUTTERITES 126251 24121 21442 1 -9
HUTTERITES 26612 5241 5452 2 -9
HUTTERITES 111161 24111 25482 1 -9
...

TPED FILE: 
22 rs2379981 0 15410792 2 2 2 2 2 2
22 rs5746647 0 15437138 2 2 2 2 2 2
...

Created on June 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import csv, numpy as np, os
from impute.data.factory import GenotypeFactory
from impute.data.io_pedigree import _FAMILY_ID, read_sample_id
from impute.tools import recode
from impute.data.constants import MISSING
from impute.tools.recode import CGI_LETTER_TO_ALLELE
# from impute.impute_test_util import HUTT_PED, GENOTYPE_SAMPLE

####################################################################################
class __FileReader(object):
    '''An abstraction of file readers.'''
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def read(self, clazz, **kwargs):
        '''Read a data set of type clazz from files/streams.'''
        raise 'Sub-classes must implement this method' 

####################################################################################
class __PlinkFileReader(__FileReader):
    '''Read a Genotype object from a PLINK file set.'''
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def read(self, clazz, **kwargs):
        '''Load genotype data. If prefix is specified, will use prefix.tfam, prefix.tped
        input file names, unless tfam and/or tped are specified (with or without the
        prefix argument), in which case they override the prefix-based names.'''

        # Read input arguments        
        prefix = kwargs.get('prefix', None)
        load_ids = kwargs.get('load_ids', True)
        tped = kwargs.get('tped', None if prefix is None else (prefix + '.tped'))
        if tped is None:
            raise ValueError('Must specify plink file prefix and/or tped file name')
        if load_ids: 
            tfam = kwargs.get('tfam', None if prefix is None else (prefix + '.tfam'))
            if tfam is None:
                raise ValueError('If loading IDs, must specify plink file prefix and/or tfam file name')
        # lazily-load data or not fetch all of it 
        lazy_load = kwargs.get('lazy_load', False)
        
        # Read TPED file in two sweeps.        
        # See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map
        
        # Read the first line in the file to determine the number of samples
        with open(tped, 'r') as f:
            reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
            line = reader.next()
            if line[-1] == '': line = line[:-1]  # Trim last item in field list of this line if it is blank
            num_items = len(line)
        
        # Read SNP metadata into a record array
        snp_dtype = [('chrom', np.uint8),  # Chromosome # containing the SNP
                     ('name', np.chararray),  # SNP name (e.g., 'rs...')
                     ('dist_cm', np.float),  # Genetic position [CENTI-Morgans!!]
                     ('base_pair', np.uint)  # Base pair position on chromosome
                    ]
        snp = np.loadtxt(tped, usecols=range(4), dtype=snp_dtype)
        # Fix the special case of a single row, where loadtxt is buggy
        if snp.size == 1: snp = np.array([tuple(snp[key] for key, _ in snp_dtype)], dtype=snp_dtype)
        
        # Read Genotype data
        if lazy_load:
            # Only pass pointer to file, to be read into a data structure that supports lazy loading
            data = tped
        else:
            # Read Genotype data into array
            data = np.genfromtxt(tped, usecols=range(4, num_items), dtype=np.byte)
            if np.size(snp) == 1: data = data.reshape([1, data.shape[0] / 2, 2])
            else: data = data.reshape([data.shape[0], data.shape[1] / 2, 2])
                
        # Load TFAM data, use only study IDs
        sample_id = np.genfromtxt(tfam, dtype=np.int)[:, 1] if load_ids else None

        # Construct object
        return GenotypeFactory.new_instance(clazz, data, snp, sample_id, lazy_load=lazy_load)

####################################################################################
class FileWriter(object):
    '''An abstraction of genotype file writers.'''
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def write(self, genotype, out, **kwargs):
        '''Write a data set of type clazz to the output stream out.'''
        raise ValueError('Sub-classes must implement this method')

####################################################################################
class __PlinkFileWriter(FileWriter):
    '''Write a Genotype object to a PLINK TPED file.'''
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def write(self, genotype, out, **kwargs):
        '''Write genotype data to the stream out in tped format. If snps is specified, 
        only those snp indices are written, otherwise all snps are written.''' 

        # Read optional arguments
        snps = kwargs.get('snps', genotype.snp_range)
        
        data, fmt = genotype.data, '%d'
        if kwargs.get('recode_cgi', False): data, fmt = recode.recode_cgi(data), '%s'
            
        # Convert SNP data to string
        snp_metadata = genotype.snp  # [np.array(['chrom', 'name', 'dist_cm', 'base_pair'])]
        snp_as_string = np.array([str(y) for x in snp_metadata for y in x]).reshape((genotype.num_snps, 4))
        for snp in snps:
            np.savetxt(out, snp_as_string[snp], fmt='%s', newline=' ')
            np.savetxt(out, data[snp, :, :], fmt=fmt, newline=' ')
            out.write('\n')
            
        if kwargs.get('sample_id_out', None):
            np.savetxt(kwargs.get('sample_id_out'), genotype.sample_id,
                       fmt=_FAMILY_ID + ' %s 0 0 1 -9 0', newline='\n')

####################################################################################
class __GaixinFileWriter(FileWriter):
    '''Write a Genotype object to file in the format requested by Gaixin for her analyses.
    This generates a text file with four columns: SNP ID, sample ID, allele 1, allele2 for all
    samples and individuals in the Genotype object.'''
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def write(self, g, out, **kwargs):
        '''Write genotype data to the stream out in Gaixin format..''' 
        np.savetxt(out, np.concatenate((np.tile(g.snp['name'], (g.num_samples, 1))
                                        .transpose().flatten()[np.newaxis].transpose(),
                                        np.tile(kwargs.get('sample_id', np.arange(g.num_samples)),
                                                (1, g.num_snps)).transpose(),
                                        g.data.reshape((g.data.shape[0] * g.data.shape[1], 2))), axis=1),
                   fmt='%s %d %d %d')

####################################################################################
class __Impute2FileWriter(FileWriter):
    '''Write a Genotype object to file in IMPUTE2 haplotype format.'''
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def write(self, genotype, out, **kwargs):
        '''Write genotype data to the stream out in IMPUTE2 haplotype format. If snps is specified, 
        only those snp indices are written, otherwise all snps are written.''' 

        # Read optional arguments; truncate data accordingly
        snps = kwargs.get('snps', genotype.snp_range)
        samples = kwargs.get('samples', np.arange(genotype.num_samples))
        flip_alleles = kwargs.get('flip_alleles', np.zeros(len(snps), dtype=bool))
        
        # Zero out partially-called genotypes (an IMPUTE2 requirement)
        data = genotype.data[snps][:, samples, :].copy()
        recode.clear_partial_calls(data)
        
        # Convert SNP data to string
        snp_metadata = genotype.snp[snps]  # [np.array(['chrom', 'name', 'dist_cm', 'base_pair'])]
        snp_as_string = np.array([str(y) for x in snp_metadata for y in x]).reshape((len(snps), 4))

        # Recode alleles + missing data
        recoding = { MISSING: '?', 1: '0', 2: '1' }
        r = np.vectorize(recoding.__getitem__)
        data_str = r(data)
        # Flip alleles
        recoding_flipped = { MISSING: '?', 1: '1', 2: '0' }
        r = np.vectorize(recoding_flipped.__getitem__)
        data_str[flip_alleles] = r(data[flip_alleles]) 
        
        for i in xrange(len(snps)):
            np.savetxt(out, snp_as_string[i], fmt='%s', newline=' ')
            np.savetxt(out, data_str[i, :, :], fmt='%s', newline=' ')
            out.write('\n')

#---------------------------------------------
# Facade Methods
#---------------------------------------------

_plink_reader = __PlinkFileReader()

_plink_writer = __PlinkFileWriter()
_gaixin_writer = __GaixinFileWriter()
_impute2_writer = __Impute2FileWriter()

def read(input_type, clazz, **kwargs):
    '''Read a Genotype object of class clazz ('genotype'/'haplotype'/'problem')
    from file of the format 'input_type'.
    
    Supported formats: input_type='plink' (PLINK format); 'npz' (our NPZ format).'''
    if input_type == 'npz':
        data = np.load(kwargs.get('file'))
        g = GenotypeFactory.new_instance(clazz, data['data'], data['snp'], sample_id=data['sample_id'])
        # If there exists a genetic map, load it. If not, don't. For backward-compatibility with older
        # Genotype npz files that didn't have the map yet
        if 'map' in data.files: g.map = data['map']
        if 'poo_phase' in data.files: g.poo_phase = data['poo_phase']
        return g
    elif input_type == 'plink':
        return _plink_reader.read(clazz, **kwargs)
    else:
        raise ValueError('Unsupported genotype input type %s' % (input,))

def write(output_type, genotype, out, **kwargs):
    '''Write genotype data to the stream out in tped format. If snps is specified, 
    only those snp indices are written, otherwise all snps are written. 
    Currently only output_type='plink' is supported: PLINK format.'''
    if output_type == 'npz': np.savez(out, data=genotype.data, snp=genotype.snp, sample_id=genotype.sample_id, map=genotype.map, poo_phase=genotype.poo_phase)
    elif output_type == 'plink': _plink_writer.write(genotype, out, **kwargs)
    elif output_type == 'gaixin': _gaixin_writer.write(genotype, out, **kwargs)
    elif output_type == 'impute2': _impute2_writer.write(genotype, out, **kwargs)
    else: raise ValueError('Unsupported genotype output type %s' % (input,))

def read_tabix(file_name, genotyped_id_file=os.environ['OBER_DATA'] + '/hutt/hutt.3chipoverlap.clean.fam'):
    '''Read a Haplotype object from an ITABIX CGI-imputed file.
    Line format: tab-delimited
    7849538    chr11    1909005    1909006    snp    T    C    dbsnp.107:rs3817198    <genotypes>
    '''
    # Load entire file into memory. It must fit, if we are to load it into a Genotype object
    d = np.loadtxt(file_name, str) 
        
    # Read SNP metadata into a record array
    snp_dtype = [('chrom', np.uint8),  # Chromosome # containing the SNP
                 ('name', np.chararray),  # SNP name (e.g., 'rs...')
                 ('dist_cm', np.float),  # Genetic position [CENTI-Morgans!!]
                 ('base_pair', np.uint)  # Base pair position on chromosome
                ]
    snp = np.array([(int(line[1][3:]), line[7], 0, int(line[3])) for line in d], dtype=snp_dtype)
    data = np.array([[(CGI_LETTER_TO_ALLELE[x[1]], CGI_LETTER_TO_ALLELE[x[2]]) for x in line[8:]] for line in d])
    hap_type = np.array([[int(x[0]) for x in line[8:]] for line in d])
    sample_id = read_sample_id(genotyped_id_file)
    # Construct object
    return GenotypeFactory.new_instance('haplotype', data, snp, sample_id, hap_type=hap_type)
