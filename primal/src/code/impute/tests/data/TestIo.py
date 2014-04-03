'''
============================================================
Test loading and saving Problem objects in PLINK and NPZ
formats.

Created on August 5, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

import tempfile, unittest, numpy as np, networkx as nx, os, traceback
from impute import impute_test_util as itu
from numpy.ma.testutils import assert_equal
from impute.phasing.phase_trivial import trivial_phaser
from impute.data import io, io_pedigree
from impute.phasing.phase_core import new_phaser_chain

class TestIo(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Load test data and expected results.'''
        unittest.TestCase.setUp(self)
        # Test phaser, doesn't matter which here
        self.phaser = new_phaser_chain([trivial_phaser()])
        # Generate a non-trivial Problem
        self.problem = itu.Templates.problem_family(itu.FAMILY13, haplotype=True)

    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_save_load_problem_npz(self):
        '''Check that saving and loading a problem from file preserves the original object.'''
        self.__test_save_and_load_problem(self.problem, 'npz')
    
    def test_save_load_problem_plink(self):
        '''Check that saving and loading a problem from file preserves the original object.'''
        self.__test_save_and_load_problem(self.problem, 'plink')

    def test_save_load_problem_with_haplotype(self):
        '''Check that saving and loading a problem from file preserves the original object.'''
        p = itu.Templates.problem_family(itu.FAMILY13, haplotype=True)
        self.phaser.run(p)
        self.__test_save_and_load_problem(p, 'npz')

    def test_save_load_problem_with_problem_info(self):
        '''Check that saving and loading a problem from file preserves the original object.'''
        self.__test_save_and_load_problem(io.read_npz(itu.FAMILY13_STAGE3), 'npz')
        self.__test_save_and_load_problem(io.read_npz(itu.FAMILY13_STAGE3), 'plink')

    def test_save_load_subproblem(self):
        '''Check saving and loading a sub-problem with selected samples.'''
        p = io.read_npz(itu.FAMILY13_STAGE3)
        pp = p.sub_problem(p.families()[0].member_list) 
        self.__test_save_and_load_problem(pp, 'npz')
        self.__test_save_and_load_problem(pp, 'plink')

    def test_save_load_pedigree_graph(self):
        '''Check that saving and loading a DiGraph from file preserves the original object.'''
        self.__test_save_and_load_graph_npz(nx.path_graph(10, nx.DiGraph()))
        self.__test_save_and_load_graph_npz(self.problem.pedigree.graph)

    def test_save_load_pedigree_plink(self):
        '''Check that saving and loading a pedigree object from file preserves the original object.'''
        p = itu.Templates.pedigree_hut()
        out_file = tempfile.TemporaryFile()
        io_pedigree.write(p, out_file)
        out_file.seek(0)
        p2 = io_pedigree.read( out_file, genotyped_id_file=itu.GENOTYPE_SAMPLE+'.tfam')
        out_file.close()
        assert_equal(p, p2, 'Saving and loading did not restore the original pedigree')

    def test_plink_to_npz(self):
        '''Check that converting plink to npz works.'''
        self.__test_convert(self.__plink_to_npz)

    def test_npz_to_plink(self):
        '''Check that converting npz to plink works.'''
        self.__test_convert(self.__npz_to_plink)

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __test_save_and_load_problem(self, problem, fmt):
        '''Check that saving and loading a problem from file preserves the original
        object.'''
        if fmt == 'npz':
            p = self.__save_and_load_problem_npz(problem)
        else:
            p = self.__save_and_load_problem_plink(problem)
        assert_equal(p.pedigree.num_genotyped, problem.pedigree.num_genotyped, 'Saving and loading did not restore problem pedigree # genotyped samples')
        assert_equal(p.pedigree.parents(1), problem.pedigree.parents(1), 'Saving and loading did not restore problem pedigree graph edge type attributes') 
        assert_equal(p.error, problem.error, 'Saving and loading did not restore the error array')
        assert_equal(p, problem, 'Saving and loading did not restore the original object')
        #out_file.close() # Would be needed with a normal file descriptor
    
    def __save_and_load_problem_npz(self, problem):
        '''Save and load a problem from NPZ file.'''
        out_file = tempfile.TemporaryFile()
        io.write_npz(problem, out_file)
        # Only needed here to simulate closing & reopening file; you will need to call 
        # out_file.close() on to prevent file locking in Windows
        out_file.seek(0)
        return io.read_npz(out_file)

    def __save_and_load_problem_plink(self, problem):
        '''Save and load a problem from PLINK file set.'''
        try:
            # Get a temporary file name
            f = tempfile.NamedTemporaryFile(delete=False)
            file_name = f.name
            f.close()
            io.write_plink(problem, file_name)
            return io.read_plink(prefix=file_name)
        finally:
            # Delete test files
            for ext in ['', '.pdg.tfam', '.tfam', '.tped', '.hap.tped', '.info']:
                os.remove(file_name + ext)

    def __test_save_and_load_graph_npz(self, x):
        '''Test save and load a Networkx DiGraph in npz format with np-array wrapping.'''
        out_file = tempfile.TemporaryFile()
        np.savez(out_file, x=np.array([nx.to_scipy_sparse_matrix(x)]))
        out_file.seek(0) # Only needed here to simulate closing & reopening file
        x2 = np.load(out_file)
        y = nx.from_scipy_sparse_matrix(x2['x'][0], nx.DiGraph())
        assert_equal(x.nodes(), y.nodes(), 'Saving and loading did not restore the original object')
        assert_equal(x.edges(), y.edges(), 'Saving and loading did not restore the original object')
        #out_file.close() # Would be needed with a normal file descriptor
    
    def __test_convert(self, convert):
        '''Check that converting from one format to another works.'''
        p = self.problem
        try:
            # Get a temporary file name
            f = tempfile.NamedTemporaryFile(delete=False)
            file_name = f.name
            f.close()
            p2 = convert(p, file_name)
            assert_equal(p.pedigree.graph.nodes(), p2.pedigree.graph.nodes(), 'Saving and loading did not restore the original Problem - in pedigree graph nodes')
            assert_equal(p.pedigree.graph.edges(), p2.pedigree.graph.edges(), 'Saving and loading did not restore the original Problem - in pedigree graph edges')
            #assert_equal(p.pedigree.sex, p2.pedigree.sex, 'Saving and loading did not restore the original Problem - in pedigree sex')
            assert_equal(p.pedigree, p2.pedigree, 'Saving and loading did not restore the original Problem - in pedigree')
            assert_equal(p.genotype, p2.genotype, 'Saving and loading did not restore the original Problem - in genotype')
            assert_equal(p.haplotype, p2.haplotype, 'Saving and loading did not restore the original Problem - in haplotype')
            assert_equal(p.info, p2.info, 'Saving and loading did not restore the original Problem - in ProblemInfo')
            assert_equal(p, p2, 'Saving and loading did not restore the original Problem')
        except Exception, e:
            traceback.print_exc()
            raise e
        finally:
            # Delete test files
            try:
                for ext in ['', '.pdg.tfam', '.tfam', '.tped', '.hap.tped', '.info', '.npz']:
                    os.remove(file_name + ext)
            except:
                pass

    def __plink_to_npz(self, p, file_name):
        '''Convert p from plink to npz format using the file set specified by file_name.'''
        npz = file_name+'.npz'
        # Save test problem in plink format
        io.write_plink(p, prefix=file_name)
        # Convert plink -> npz
        io.plink_to_npz(file_name, npz)
        # Load npz and check that the problem object didn't change
        p2 = io.read_npz(npz)
        return p2

    def __npz_to_plink(self, p, file_name):
        '''Convert p from npz to plink format using the file set specified by file_name.'''
        npz = file_name+'.npz'
        # Save test problem in plink format
        io.write_npz(p, npz)
        # Convert plink -> npz
        io.npz_to_plink(npz, file_name)
        # Load npz and check that the problem object didn't change
        p2 = io.read_plink(prefix=file_name)
        return p2
