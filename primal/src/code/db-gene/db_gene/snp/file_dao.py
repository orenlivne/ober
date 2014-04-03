'''
============================================================
DAO implementation that reads from text files, not a
database.

Created on July 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, csv, linecache, itertools as it
from scipy.sparse.csr import csr_matrix
from scipy.sparse import triu

#---------------------------------------------
# Methods
#---------------------------------------------
def standardize_kinship_format(in_file, out_file):
    '''Read the upper-triangular part of the kinship matrix from in_file. Convert to the
    format expected by KinshipDao and write to out_file.
    
    Note: loads the entire file into memory.'''
    data = np.loadtxt(in_file, usecols=[2])
    n = int(((8 * len(data) + 1) ** 0.5 - 1) / 2)
    idx = np.array(list(it.chain.from_iterable(xrange(k, n) for k in xrange(n))))
    idx_ptr = np.concatenate(([0], np.cumsum(xrange(n, 0, -1))))
    A = csr_matrix((np.maximum(data, 1e-16), idx, idx_ptr), shape=(n, n))
    with open(out_file, 'wb') as f:
        f.write(' '.join(it.islice((x[1] for x in csv.reader(open(in_file, 'rb'), delimiter='\t')), n)) + '\n')
    with open(out_file, 'ab') as f:
        np.savetxt(f, (A + triu(A, 1).transpose()).data, fmt='%.16f')
    return A

####################################################################################
class ChromDao(object):
    '''Retrieves chromosome metadata.''' 
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    # Chromosome size [base pair]
    TOTAL_BP = [247199719L, 242751149L, 199446827L, 191263063L, 180837866L, 170896993L,
                158821424L, 146274826L, 140442298L, 135374737L, 134452384L, 132289534L,
                114127980L, 106360585L, 100338915L, 88822254L, 78654742L, 76117153L,
                63806651L, 62435965L, 46944323L, 49528953L]

    # Chromosome size [base pair] (the part typed by the affy chip)
    TOTAL_BP_TYPED = [248362901L, 242642248L, 197726641L, 190836115L, 180580121L, 170740586L,
                      159059338L, 146102846L, 140945683L, 135317340L, 134731703L, 133571132L,
                      96732839L, 87133853L, 82183303L, 90114806L, 80939731L, 77853695L, 60898154L,
                      62837106L, 33321451L, 34672141L]

    # Chromosome size [cM]
    TOTAL_CM = [286.279234,
                268.839622,
                223.361095,
                214.688476,
                204.089357,
                192.039918,
                187.220500,
                168.003442,
                166.359329,
                181.144008,
                158.218650,
                174.679023,
                125.706316,
                120.202583,
                141.860238,
                134.037726,
                128.490529,
                117.708923,
                107.733846,
                108.266934,
                62.786478,
                74.109562]
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def total_bp(self):
        '''Return the total # base pairs on each chromosome. Ordered by chromosome #.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        return ChromDao.TOTAL_BP

    def total_bp_typed(self):
        '''Return the total # base pairs on each chromosome in the region typed by the 
        Hutterites Affymetrix chip. Ordered by chromosome #.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        return ChromDao.TOTAL_BP_TYPED

    def total_cm(self):
        '''Return the total # base pairs on each chromosome. Ordered by chromosome #.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        return ChromDao.TOTAL_CM

####################################################################################
class IdCoefDao(object):
    '''Retrieves Hutterites identity coefficients from a file. The file is assumed to contain
    n^2 entries.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, file_name, delimiter=' '):
        '''Initialize from a data file file_name.'''
        with open(file_name, 'r') as f:
            reader = csv.reader(f, delimiter=delimiter, skipinitialspace=True)
            self.index = {}
            prev_id1 = None
            # Loop over lines until the first column changes
            for i, line in enumerate(reader):
                id1, id2 = int(line[0]), int(line[1])
                if prev_id1 and prev_id1 != id1: break
                prev_id1 = id1
                self.index[id2] = i
            self.file_name = file_name
            self.n = len(self.index)
            self.delimiter = delimiter
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of kinship pairs (n^2 where n=#samples).'''
        return self.n ** 2
        
    def id_coefs(self, id1, id2):
        '''Load the condensed identity coefficients (lambda, (Delta1,...,Delta9)) between id1 and id2.'''
        line = linecache.getline(self.file_name, self.__line_number(id1, id2))
        items = line.rstrip('\n').split(self.delimiter)
        return float(items[2]), np.array([float(items[x]) for x in xrange(3, 12)])
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __line_number(self, id1, id2):
        '''Return the 1-based line number of the pair (id1,id2). This is look-up in a full matrix.'''
        return self.n * self.index[id1] + self.index[id2] + 1

####################################################################################
class KinshipDao(object):
    '''Retrieves Hutterites kinship coefficients from a file. The file is assumed to contain
    n^2 entries.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, file_name):
        '''Initialize from a data file file_name.'''
        self.file_name = file_name
        # Read sample ID index from first line in the file
        with open(file_name, 'r') as f:
            line = f.next()
            ids = np.fromstring(line, dtype=int, sep=' ')
            self.index = dict((v, k) for k, v in enumerate(ids))
            self.n = len(self.index)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of kinship pairs (n^2 where n=#samples).'''
        return self.n ** 2
        
    def kinship(self, id1, id2):
        '''Load the kinship coefficient between id1 and id2.'''
        line = linecache.getline(self.file_name, self.__line_number(id1, id2)) 
        return float(line.rstrip('\n'))
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __line_number(self, id1, id2):
        '''Return the 1-based line number of the pair (id1,id2), plus one more line for the
        index header. This is look-up in a full matrix.'''
        return self.n * self.index[id1] + self.index[id2] + 2

####################################################################################
class Daos(object):
    '''DAO mother object.'''
    def __init__(self, **kwargs):
        self.configure(**kwargs)

    def configure(self, **kwargs):
        # Synchronize access to url
        self.chrom_dao = ChromDao()

# Global access to DAO, default configuration
DEFAULT_FILE_DAOS = Daos()
