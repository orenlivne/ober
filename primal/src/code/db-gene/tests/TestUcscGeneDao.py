'''
============================================================
Test retrieving Gene information from a local mirror of the
UCSC gene browser. 

Created on November 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest
from numpy.ma.testutils import assert_equal
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.sql.expression import and_
from db_gene.ucsc.entities import Gene
from db_gene.op.functions import greatest, least
from db_gene.ucsc.ucsc_dao import GeneDao
from db_gene import DEFAULT_URL

class TestUcscGeneDao(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Create an in-memory database.'''
        self.engine = create_engine(DEFAULT_URL)
        self.dao = GeneDao(DEFAULT_URL)
        # Base.metadata.create_all(self.engine) 
        
    def tearDown(self):
        '''Drop the database.'''
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_num_records(self):
        '''Test getting the total # of refGene records.'''
        # As of 13-NOV-2012
        self.assertTrue(self.dao.num_records() >= 43766, 'Wrong row count from refGene table')

    def test_intersecting_genes(self):
        '''Test retrieving genes intersecting a base-pair position segment.'''
        # print repr(User.__table__) #@UndefinedVariable
        Session = sessionmaker(bind=self.engine)
        session = Session()
        
        # I believe all of the following are equivalent statements and all are prepared statements
        chrom = 1
        (start, end) = (62151773, 62516683)
        expected = [('NM_032027', 'TM2D1'), ('NM_176877', 'INADL')]
        
        assert_equal([(gene.name, gene.name2) for gene in session.query(Gene).
                      filter(and_(Gene.chrom == 'chr%s' % (chrom,),
                                  greatest(Gene.txStart, start) <= 
                                  least(Gene.txEnd, end))).
                      order_by(Gene.txStart)], expected)

        assert_equal([(gene.name, gene.name2) for gene in session.query(Gene).
                      filter('chrom = :chrom and greatest(txStart, :start) <= least(txEnd, :end)').
                      params(chrom='chr%s' % (chrom,), start=start, end=end).
                      order_by(Gene.txStart)], expected)

        assert_equal([(gene.name, gene.name2) for gene in 
                      self.dao.intersecting_genes(chrom, start, end)], expected)
        
        session.close()
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
