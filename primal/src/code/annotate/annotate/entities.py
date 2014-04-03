'''
============================================================
Hutterites annotation database entities.

Created on April 4, 2014
@author: Oren Livne <livne@uchicago.edu>
@see http://docs.sqlalchemy.org/en/latest/orm/tutorial.html
============================================================
'''
from sqlalchemy import Column, Integer, String, Sequence, BigInteger, Float, Boolean
from sqlalchemy.ext.declarative import declarative_base

# Common base class of all entities 
Base = declarative_base()

####################################################################################
class Variant(Base):
    '''Variant annotation.
    See http://workshops.arl.arizona.edu/sql1/sql_workshop/mysql/ucscdatabase.html#querying-the-refgene-table
    '''
    __tablename__ = 'hutt'

    #id = Column(Integer, Sequence('variant_id_seq'), primary_key=True, name='variant_id')  # @ReservedAssignment
    record_id = Column(BigInteger, primary_key=True, unique=True)
    maf_imputed = Column(Float) # Minor allele frequency in the imputed Hutterites
    is_qc = Column(Boolean) # Did variant pass QC?
    var_region = Column(String(255))  # Annotation - exonic/intronic/...
    var_mutation = Column(String(255))  # Annotation - missense/nonsense/...

    def __repr__(self):
        return "Gene(id=%d, '%s')" % (self.id, self.name2,)
