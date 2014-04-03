'''
============================================================
UCSC browser database entities.

Created on October 30, 2012
@author: Oren Livne <livne@uchicago.edu>
@see http://docs.sqlalchemy.org/en/latest/orm/tutorial.html
============================================================
'''
from sqlalchemy import Column, Integer, String, Sequence, BigInteger  # , ForeignKey
from sqlalchemy.ext.declarative import declarative_base
# from sqlalchemy.orm import relationship, backref #@UnresolvedImport

# Common base class of all entities 
Base = declarative_base()

####################################################################################
class Gene(Base):
    '''Basic gene annotation (e.g., start and stop position).
    See http://workshops.arl.arizona.edu/sql1/sql_workshop/mysql/ucscdatabase.html#querying-the-refgene-table
    '''
    __tablename__ = 'refGene'

    id = Column(Integer, Sequence('refGene_id_seq'), primary_key=True, name='refGene_id')  # @ReservedAssignment
    name = Column(String(255))
    chrom = Column(Integer)
    txStart = Column(Integer)
    txEnd = Column(Integer)
    name2 = Column(String(255))  # Altername name

#    def __init__(self, name, fullname, password):
#        self.name = name
#        self.fullname = fullname
#        self.password = password

    def __repr__(self):
        return "Gene(id=%d, '%s')" % (self.id, self.name2,)

####################################################################################
class Snp(Base):
    '''All-SNP annotations (e.g., rs-number, base pair position) - .'''
    DBSNP_VERSION = '137'  # Latest available version; using hg19 = build37
    __tablename__ = 'snp%s' % (DBSNP_VERSION,)

    id = Column(Integer, Sequence('snp%s_id_seq' % (DBSNP_VERSION,)), primary_key=True, name='snp%s_id' % (DBSNP_VERSION,))  # @ReservedAssignment
    name = Column(String, primary_key=True, unique=True)
    chrom = Column(Integer, index=True)
    bp = Column(BigInteger, name='chromEnd')

#    def __init__(self, chrom, name, bp):
#        '''Populate all data fields.'''
#        self.chrom = chrom
#        self.name = name
#        self.bp = int(bp)

    def __repr__(self):
        return "Snp[%s:'%s', %d]" % (self.chrom, self.name, self.bp)
