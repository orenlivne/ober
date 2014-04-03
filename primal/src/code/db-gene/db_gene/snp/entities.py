'''
============================================================
Our Genome db (db gene) database entities.

Created on November 8, 2012
@author: Oren Livne <livne@uchicago.edu>
@see http://docs.sqlalchemy.org/en/latest/orm/tutorial.html
============================================================
'''
from sqlalchemy import Column, Integer, SmallInteger, BigInteger, String, Sequence, Float  # , ForeignKey
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.schema import Index
# from sqlalchemy.orm import relationship, backref #@UnresolvedImport

# Common base class of all entities 
Base = declarative_base()

####################################################################################
class Chromosome(Base):
    '''Chromosome-level annotations (e.g., start and stop positions).'''
    __tablename__ = 'chromosome'

    id = Column(Integer, Sequence('chromosome_id_seq'), primary_key=True, name='id')  # @ReservedAssignment
    number = Column(SmallInteger)  # numeric coding: X=23, Y=24
    name = Column(String(2))  # Chromosome name (1-22|X|Y)
    num_genes = Column(Integer)
    total_bp = Column(BigInteger)
    sequenced_bp = Column(BigInteger)

    def __init__(self, number, name, num_genes, total_bp, sequenced_bp):
        '''Populate all data fields.'''
        self.number = int(number)
        self.name = name
        self.num_genes = int(num_genes)
        self.total_bp = long(total_bp)
        self.sequenced_bp = long(sequenced_bp)

    def __repr__(self):
        return 'Chromosome[%s]' % (self.name,)

####################################################################################
class Snp(Base):
    '''SNP annotations (e.g., rs-number, base pair position) only for SNPs of interest (e.g., the SNPs
    for which the Hutterites are Affy-genotyped).'''
    __tablename__ = 'snp'

    id = Column(Integer, Sequence('snp_id_seq'), primary_key=True, name='id')  # @ReservedAssignment
    chrom = Column(SmallInteger, index=True, nullable=False)
    name = Column(String(30), index=True, unique=True, nullable=False)
    bp = Column(BigInteger, index=True, nullable=False)  # Base-pair position
    genetic_pos = Column(Float)  # Genetic position

    def __init__(self, chrom, name, bp):
        '''Populate all data fields.'''
        self.chrom = chrom
        self.name = name
        self.bp = int(bp)

    def __repr__(self):
        return 'Snp[%s, %d]' % (self.name, self.bp)

####################################################################################
class Ld(Base):
    '''LD measure between two Snp entities.'''
    __tablename__ = 'ld'

    id = Column(Integer, Sequence('ld_id_seq'), primary_key=True, name='id')  # @ReservedAssignment
    chrom = Column(SmallInteger, nullable=False)
    snp1 = Column(String(30), nullable=False)
    snp2 = Column(String(30), nullable=False)
    r2 = Column(Float, nullable=False)  # r^2 measure

    def __init__(self, chrom, snp1, snp2, r2):
        '''Populate all data fields.'''
        self.chrom = chrom
        self.snp1 = snp1
        self.snp2 = snp2
        self.r2 = float(r2)

    def __repr__(self):
        return 'Ld[%s-%s: %.2f]' % (self.snp1, self.snp2, self.r2)

####################################################################################
class Kinship(Base):
    '''Condensed identity coefficients between each two Hutterite samples.'''
    __tablename__ = 'kinship'
    # Define indexes
    __table_args__ = (Index('index_id1_id2', 'id1', 'id2'),)
 
    id = Column(Integer, Sequence('kinship_id_seq'), primary_key=True, name='id')  # @ReservedAssignment
    id1 = Column(Integer, nullable=False)
    id2 = Column(Integer, nullable=False)
    lam = Column(Float, nullable=False)
    delta1 = Column(Float, nullable=False)
    delta2 = Column(Float, nullable=False)
    delta3 = Column(Float, nullable=False)
    delta4 = Column(Float, nullable=False)
    delta5 = Column(Float, nullable=False)
    delta6 = Column(Float, nullable=False)
    delta7 = Column(Float, nullable=False)
    delta8 = Column(Float, nullable=False)
    delta9 = Column(Float, nullable=False)

    def __init__(self, id1, id2, lam, delta):
        '''Populate all data fields.'''
        self.id1 = id1
        self.id2 = id2
        self.lam = float(lam)
        # To-do: replace by one field, self.delta
        self.delta1 = float(delta[0])
        self.delta2 = float(delta[1])
        self.delta3 = float(delta[2])
        self.delta4 = float(delta[3])
        self.delta5 = float(delta[4])
        self.delta6 = float(delta[5])
        self.delta7 = float(delta[6])
        self.delta8 = float(delta[7])
        self.delta9 = float(delta[8])

    def __repr__(self):
        return 'Kinship[%s-%s: lambda=%.2f, Delta=%s]' % (self.id1, self.id2, self.lam, repr(self.delta))

    @property
    def delta(self):
        return [self.delta1, self.delta2, self.delta3, self.delta4, self.delta5, self.delta6, self.delta7, self.delta8, self.delta9]
