'''
============================================================
An immutable object representing a haploid organism (e.g.,
a person). Serves as the node type in a Pedigree's graph.

Created on August 15, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from impute.data.constants import MISSING

class Person(object):
    '''An immutable object representing a haploid organism (e.g., a person). Serves as the
    node type of a Pedigree graph.'''
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    # Affection status missing data code for quantitative traits. Matches the plink standard
    # http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
    PHENOTYPE_MISSING = -9

    class SEX:
        UNKNOWN, MALE, FEMALE = range(3)
    class TYPE:
        NOT_GENOTYPED, GENOTYPED, DUMMY, MARRIAGE = range(4)

    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, sample_id, sex=SEX.UNKNOWN, phenotype=PHENOTYPE_MISSING, 
                 node_type=TYPE.DUMMY, father=MISSING, mother=MISSING):
        '''
        Construct a person with a unique identifier _sample_id.
        ID conventions: see documentation of function id_to_sex().
        '''
        # Public properties
        self.sample_id  = sample_id
        self.sex        = sex
        self.node_type  = node_type
        self.phenotype  = phenotype
        self.father     = father
        self.mother     = mother
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Person[id=%d, sex=%d, phenotype=%d, type=%d, father=%d, mother=%d]' % \
            (self.sample_id, self.sex, self.phenotype, self.node_type, self.father, self.mother) 

    def __key(self):
        '''Hash key.'''
        return self.sample_id

    def __eq__(self, other):
        '''Equality of objects.'''
        return self.__key() == other.__key()

    def __ne__(self, other):
        '''Inequality of objects.'''
        return self.__key() != other.__key()

    def __cmp__(self, other):
        return cmp(self.sample_id, other.sample_id)
    
    def copy(self):
        '''Return a deep copy of this object.'''
        return Person(sample_id=self.sample_id, sex=self.sex, phenotype=self.phenotype,
                      node_type=self.node_type, father=self.father, mother=self.mother) 
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------