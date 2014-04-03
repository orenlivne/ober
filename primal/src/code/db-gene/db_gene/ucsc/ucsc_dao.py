'''
============================================================
UCSC Genome browser - gene annotation Data Access Object (DAO).

Created on November 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from db_gene.dao import Dao # Must be the first import
from sqlalchemy import func
from sqlalchemy.sql.expression import and_
from db_gene.ucsc.entities import Gene, Snp
from db_gene import DEFAULT_URL
from db_gene.op.functions import greatest, least

####################################################################################
class GeneDao(Dao):
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        super(GeneDao, self).__init__(url, **kwargs)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of records in the refGene database.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        session = self.Session()
        result = self.Session().query(func.count(Gene.id)).scalar()
        session.close()
        return result

    def intersecting_genes(self, chrom, start, end):
        '''Return a collection of Gene entities that overlap the base-pair [start,end] on chromosome
        chrom.'''
        # print 'intersecting_genes', chrom, start, end
        session = self.Session()
        result = list(session.query(Gene).filter(and_(Gene.chrom == 'chr%s' % (chrom,),
                                                      greatest(Gene.txStart, start) <= least(Gene.txEnd, end))).\
                                                      order_by(Gene.txStart))
        session.close()
        return result

####################################################################################
class SnpDao(Dao):
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        super(SnpDao, self).__init__(url, **kwargs)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of records in the refGene database.'''
        # return self.engine.execute("select count(*) from refGene").scalar()
        session = self.Session()
        result = self.Session().query(func.count(Snp.id)).scalar()
        session.close()
        return result

    def chrom_snps(self, chrom):
        '''Return the list of SNPs on chromosome chrom.'''
        # print 'intersecting_genes', chrom, start, end
        session = self.Session()
        result = list(session.query(Snp).filter(Snp.chrom == 'chr%s' % (chrom,)))
        session.close()
        return result

    def get_snps(self, snp_names):
        '''Return information on a list of SNPs.'''
        # print 'intersecting_genes', chrom, start, end
        session = self.Session()
        result = list(session.query(Snp).filter(Snp.name.in_(snp_names)))
        session.close()
        return result

####################################################################################
class Daos(Dao):
    '''DAO mother object.'''
    def __init__(self, url, **kwargs):
        self.configure(url, **kwargs)

    def configure(self, url, **kwargs):
        # Synchronize access to url
        self.gene_dao = GeneDao(url, **kwargs)
        self.snp_dao = SnpDao(url, **kwargs)

# Global access to DAO, default configuration. Lazily-initialized, to accommodate systems without mysql
__DEFAULT_UCSC_DAOS = None

def DEFAULT_UCSC_DAOS():
    global __DEFAULT_UCSC_DAOS
    if not __DEFAULT_UCSC_DAOS:
        __DEFAULT_UCSC_DAOS = Daos(DEFAULT_URL)
    return __DEFAULT_UCSC_DAOS
