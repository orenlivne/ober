'''
============================================================
SNP_DB Genome browser local mirror database - SNP annotation
Data Access Object (DAO).

Created on November 12, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from db_gene.dao import Dao # Must be the first import
import networkx as nx, numpy as np
from sqlalchemy import func
from sqlalchemy.sql.expression import and_, or_
from db_gene.snp.entities import Ld, Snp, Chromosome, Kinship
from db_gene import DEFAULT_URL

####################################################################################
class ChromDao(Dao):
    '''Retrieves chromosome metadata.''' 
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        super(ChromDao, self).__init__(url, **kwargs)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def total_bp(self):
        '''Return the total # base pairs on each chromosome. Ordered by chromosome #.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        session = self.Session()
        result = [instance.total_bp for instance in session.query(Chromosome).order_by(Chromosome.number)]
        session.close()
        return result

####################################################################################
class LdDao(Dao):
    '''Retrieves the SNP LD graph.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        super(LdDao, self).__init__(url, **kwargs)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_snp_records(self):
        '''Return the total number of records in the SNP table.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        session = self.Session()
        result = session.query(func.count(Snp.id)).scalar()
        session.close()
        return result

    def num_ld_records(self):
        '''Return the total number of records in the LD table.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        session = self.Session()
        result = session.query(func.count(Ld.id)).scalar()
        session.close()
        return result

    def ld_graph(self, chrom, snps=None, threshold=0.0):
        '''Load and convert the LD graph of chromosome chrom to a networkx Graph object. Only
        loads edges with r^2 >= threshold. If snps is not None, yields the sub-graph of those
        SNP names only.'''
        # print 'intersecting_genes', chrom, start, end
        session = self.Session()
        if threshold > 1e-15:
            condition = and_(Ld.chrom == '%d' % (chrom,), Ld.r2 >= threshold)
        else:
            condition = Ld.chrom == '%d' % (chrom,)
        if snps is not None:
            condition = and_(condition, Ld.snp1.in_(snps), Ld.snp2.in_(snps))
        lds = list(session.query(Ld).filter(condition))
        session.close()
        g = nx.Graph()
        g.add_weighted_edges_from((ld.snp1, ld.snp2, ld.r2) for ld in lds)
        return g
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    
####################################################################################
class IdCoefDao(Dao):
    '''Retrieves Hutterites kinship coefficients.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        super(IdCoefDao, self).__init__(url, **kwargs)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of kinship pairs.'''
        session = self.Session()
        result = session.query(func.count(Kinship.id)).scalar()
        session.close()
        return result
        
#    def all_id_coefs(self):
#        '''An iterator over the condensed identity coefficients (lambda, (Delta1,...,Delta9)) between all
#        pairs id1 and id2, sorted by id1, then by id2. 
#        Note: database only stores a non-redundant half of the pairs with id1 >= id2.'''
#        session = self.Session()
#        for result in session.query(Kinship).order_by(asc(Kinship.id1)).order_by(asc(Kinship.id2)).yield_per(100):
#            yield result.id1, result.id2, \
#            result.lam, [result.delta1, result.delta2, result.delta3,
#                         result.delta4, result.delta5, result.delta6,
#                         result.delta7, result.delta8, result.delta9]
#        session.close()

    def id_coefs(self, id1, id2):
        '''Load the condensed identity coefficients (lambda, (Delta1,...,Delta9)) between id1 and id2.'''
        # Database only stores a non-redundant half of the pairs with id1 >= id2
        #if id1 > id2:
        #    id1, id2 = id2, id1
        #print id1,id2
        session = self.Session()
        result = list(session.query(Kinship)
                      .filter(or_(and_(Kinship.id1 == '%d' % (id1,), Kinship.id2 == '%d' % (id2,)),
                                  and_(Kinship.id2 == '%d' % (id1,), Kinship.id1 == '%d' % (id2,)))))[0]
        session.close()
        return result.lam, np.array([
                result.delta1, result.delta2, result.delta3,
                result.delta4, result.delta5, result.delta6,
                result.delta7, result.delta8, result.delta9,
                ])

####################################################################################
class SnpDao(Dao):
    '''Retrieves and writes our SNP metadata.''' 
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        super(SnpDao, self).__init__(url, **kwargs)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of records in the SNP table.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        session = self.Session()
        result = session.query(func.count(Snp.id)).scalar()
        session.close()
        return result

    def delete(self, chrom):
        '''Delete all SNPS of a certain chromosome.'''
        session = self.Session()
        session.query(Snp).filter(Snp.chrom == chrom).delete()
        session.commit()
        session.close()

    def get(self, chrom):
        '''Return all SNPs on chromosome chrom.'''
        # print 'intersecting_genes', chrom, start, end
        session = self.Session()
        result = list(session.query(Snp).filter(Snp.chrom == '%d' % (chrom,)))
        session.close()
        return result

    def get_snps(self, snp_names):
        '''Return information on a list of SNPs.'''
        # print 'intersecting_genes', chrom, start, end
        session = self.Session()
        result = list(session.query(Snp).filter(Snp.name.in_(snp_names)))
        session.close()
        return result

    def get_snps_iter(self, snp_names):
        '''Return information on a list of SNPs. Iterator that pages through the snp_names set.
        Suitable for large snp_names lists'''
        # print 'intersecting_genes', chrom, start, end
        session = self.Session()
        buf_size = 100
        for i in xrange(0, len(snp_names), buf_size):
            for x in session.query(Snp).filter(Snp.name.in_(snp_names[i:i + buf_size])):
                yield x
        session.close()

    def save(self, snps):
        '''Save a list of Snp entities.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        session = self.Session()
        session.add_all(snps)
        session.commit()
        session.close()

####################################################################################
class Daos(Dao):
    '''DAO mother object.'''
    def __init__(self, url, **kwargs):
        self.configure(url, **kwargs)

    def configure(self, url, **kwargs):
        # Synchronize access to url
        self.chrom_dao = ChromDao(url, **kwargs)
        self.ld_dao = LdDao(url, **kwargs)
        self.idcoef_dao = IdCoefDao(url, **kwargs)
        self.snp_dao = SnpDao(url, **kwargs)

# Global access to DAO, default configuration. Lazily-initialized, to accommodate systems without mysql
__DEFAULT_SNP_DB_DAOS = None
def DEFAULT_SNP_DB_DAOS():
    global __DEFAULT_SNP_DB_DAOS
    if not __DEFAULT_SNP_DB_DAOS:
        __DEFAULT_SNP_DB_DAOS = Daos(DEFAULT_URL)
    return __DEFAULT_SNP_DB_DAOS
