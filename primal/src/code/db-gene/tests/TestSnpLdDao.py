'''
============================================================
Test retrieving Gene information from a local mirror of the
UCSC gene browser. 

Created on November 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, networkx as nx, numpy as np, db_gene, tempfile
from numpy.ma.testutils import assert_equal
from sqlalchemy import create_engine

class TestSnpLdDao(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Use a localhost UCSC copy.'''
        self.engine = create_engine(db_gene.DEFAULT_URL)
        self.ld_dao = db_gene.snp.snp_db_dao.DEFAULT_SNP_DB_DAOS().ld_dao
        self.snp_dao = db_gene.ucsc.ucsc_dao.DEFAULT_UCSC_DAOS().snp_dao
        self.my_snp_dao = db_gene.snp.snp_db_dao.DEFAULT_SNP_DB_DAOS().snp_dao
        # Base.metadata.create_all(self.engine) 
        
    def tearDown(self):
        '''Drop the database.'''
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_num_records(self):
        '''Test getting the total # of SNP records.'''
        assert_equal(self.ld_dao.num_ld_records(), 1394244, 'Wrong row count from SNP table')
        assert_equal(self.ld_dao.num_snp_records(), 271486, 'Wrong row count from SNP table')

    def test_ld_graph(self):
        '''Test retrieving SNPs by chromosome.'''
        chrom = 22
        g = self.ld_dao.ld_graph(chrom)
        assert_equal(g.number_of_nodes(), 2958, 'Wrong number of SNPs on chromosome %d' % (chrom,))

    def test_ld_graph_compute_frames(self):
        '''Test retrieving SNPs by chromosome and calculating frames (=independent SNP sets).'''
        chrom = 22
        g = self.ld_dao.ld_graph(chrom)
        blocks = nx.connected_components(g)
        assert_equal([len(x) for x in blocks],
                     [1172, 445, 229, 137, 63, 61, 52, 34, 30, 28, 25, 24, 21, 17, 15, 12, 12, 11, 11,
                      11, 10, 10, 9, 9, 9, 8, 8, 8, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5,
                      5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3,
                      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                      3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                      2, 2, 2, 2, 2, 2, 2],
                     'Unexpected LD block size')
        
        block_frames = list(db_gene.snp.ld_graph._compute_block_frames(g, self.snp_dao))
        assert_equal(np.array([len(x) for x in block_frames[0].itervalues()]),
                     [292, 392, 287, 100, 51, 20, 14, 6, 3, 2, 2, 1, 1, 1],
                     'Unexpected independent SNP frame sizes within LD block 0')
        
        frames = db_gene.snp.ld_graph._compute_frames(g, self.snp_dao) 
        frame_size = [len(x) for x in frames]
        assert_equal(sum(frame_size), g.number_of_nodes(), 'Wrong total number of SNPs in all frames')
        assert_equal(frame_size, [392, 292, 287, 151, 122, 105, 100,
                                  90, 89, 89, 89, 89, 89, 89, 89,
                                  89, 89, 89, 89, 88, 88, 88, 88, 88],
                     'Wrong final frame sizes')

    def test_frames(self):
        '''Test splitting a SNP list into frames (main method called by client codes).'''
        # Create test data (using ld_graph as well, but could come from a different source)
        chrom = 22
        independent_snp_names = ['rs0', 'rs1', 'zzzzzzz2'] # Are not in LD to other SNPs in snp_names
        snp_names = self.ld_dao.ld_graph(chrom).nodes()
        
        frames = db_gene.snp.ld_graph.frames(chrom, independent_snp_names + snp_names, self.snp_dao, self.ld_dao) 
        frame_size = [len(x) for x in frames]
        assert_equal(sum(frame_size), len(snp_names) + len(frames) * len(independent_snp_names),
                     'Wrong total number of SNPs in all frames')
        assert_equal(frame_size, np.array([395, 295, 290, 154, 125, 108, 103, 93, 92, 92, 92, 92, 92,
                                           92, 92, 92, 92, 92, 92, 91, 91, 91, 91, 91]),
                     'Wrong final frame sizes')

    def test_frames_with_different_dao(self):
        '''Test splitting a SNP list into frames (main method called by client codes) using
        our snp table instead of the snp135 UCSC table.'''
        # Create test data (using ld_graph as well, but could come from a different source)
        chrom = 22
        independent_snp_names = ['rs0', 'rs1', 'zzzzzzz2'] # Are not in LD to other SNPs in snp_names
        snp_names = self.ld_dao.ld_graph(chrom).nodes()
        
        frames = db_gene.snp.ld_graph.frames(chrom, independent_snp_names + snp_names, self.my_snp_dao, self.ld_dao) 
        frame_size = [len(x) for x in frames]
        assert_equal(sum(frame_size), len(snp_names) + len(frames) * len(independent_snp_names),
                     'Wrong total number of SNPs in all frames')
        assert_equal(frame_size, np.array([393, 295, 291, 154, 125, 108, 104, 93, 92, 92, 92, 92, 92,
                                           92, 92, 92, 92, 92, 92, 91, 91, 91, 91, 91]),
                     'Wrong final frame sizes')
                
    def test_frames_genotype_sample(self):
        '''Same as test_frames(), but with dat aset as in impute's TestPhasePipeline test suite.'''
        chrom = 22
        snp_names = np.array(['rs9605923', 'rs5747999', 'rs5746679', 'rs11089263', 'rs11089264',
                              'rs2845377', 'rs16984825', 'rs1654'])
        
        # Sort snp_names by anything (e.g., name) so that we can quickly search for items
        # within that list 
        orig_indices = snp_names.argsort()
        ld_g = self.ld_dao.ld_graph(chrom, snps=snp_names)
        # SNPs that are independent of all other SNPs 
        ind = list(set(snp_names) - set(ld_g.nodes()))
        independent = np.sort(orig_indices[np.searchsorted(snp_names[orig_indices], ind)])
        frames = db_gene.snp.ld_graph._compute_frames(ld_g, self.snp_dao)
        frame_index = [np.sort(np.concatenate((independent, orig_indices[np.searchsorted(snp_names[orig_indices], list(f))]))) for f in frames]
          
        assert_equal(sorted(list(reduce(set.union, frame_index, set([])))), np.arange(len(snp_names)), 'Union of frames should be the original set')
        assert_equal(frame_index, [[0, 3, 4, 5, 6, 7], [0, 1, 2, 5, 6, 7]], 'Wrong frames')
        
    def test_save_and_load_frames(self):
        '''Check that saving and loading frames from a text file preserves them.'''
        chrom = 22
        snp_names = np.array(['rs9605923', 'rs5747999', 'rs5746679', 'rs11089263', 'rs11089264',
                              'rs2845377', 'rs16984825', 'rs1654'])
        frames = db_gene.snp.ld_graph.Frames((chrom, x) for x in db_gene.snp.ld_graph.frames(chrom, snp_names, self.snp_dao, self.ld_dao))
        
        out_file = tempfile.TemporaryFile()
        db_gene.snp.ld_graph.write_frames(frames, out_file)
        # Only needed here to simulate closing & reopening file; you will need to call 
        # out_file.close() on to prevent file locking in Windows
        out_file.seek(0)
        loaded_frames = db_gene.snp.ld_graph.read_frames(out_file)
        out_file.close()
        
        assert_equal(loaded_frames, frames, 'Saving and loading did not restore original frames')
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
