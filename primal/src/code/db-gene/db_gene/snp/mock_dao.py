'''
============================================================
Mock DAO implementation. Good for unit testing.

Created on July 9, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np

####################################################################################
class ChromDao(object):
    '''Retrieves chromosome metadata.''' 
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    # Chromosome size
    TOTAL_BP = [
                247199719,
                242751149,
                199446827,
                191263063,
                180837866,
                170896993,
                158821424,
                146274826,
                140442298,
                135374737,
                134452384,
                132289534,
                114127980,
                106360585,
                100338915,
                88822254,
                78654742,
                76117153,
                63806651,
                62435965,
                46944323,
                49528953]
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def total_bp(self):
        '''Return the total # base pairs on each chromosome. Ordered by chromosome #.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        return ChromDao.TOTAL_BP

####################################################################################
class IdCoefDao(object):
    '''Retrieves Hutterites kinship coefficients from a file. The file is assumed to contain
    n^2 '''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, n):
        self.n = n
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of kinship pairs (n^2 where n=#samples).'''
        return self.n ** 2
        
    def id_coefs(self, id1, id2):
        '''Load the condensed identity coefficients (lambda, (Delta1,...,Delta9)) between id1 and id2.'''
        return (0.5, np.array([0, 0, 0, 0, 0, 0, 0.25, 0.5, 0.25]))  # Everyone is a sib

####################################################################################
class Daos(object):
    '''DAO mother object.'''
    def __init__(self, **kwargs):
        self.configure(**kwargs)

    def configure(self, **kwargs):
        # Synchronize access to url
        self.chrom_dao = ChromDao()
        # self.idcoef_dao = IdCoefDao(url, **kwargs)

# Global access to DAO, default configuration
DEFAULT_FILE_DAOS = Daos()
