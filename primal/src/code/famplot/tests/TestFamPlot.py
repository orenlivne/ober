'''
============================================================
Test pedigree plot using the famplot script (translation
of pedfiddler v0.6). 

Created on August 17, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, inspect, tempfile, unittest, test_util, numpy as np, famplot as fp
from numpy.ma.testutils import assert_equal

def abs_path(path):
    '''Path to test file.'''
    return os.path.dirname(inspect.getfile(inspect.currentframe())) + '/' + path

class TestFamPlot(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    '''pedfiddler's base example.'''
    __PEDIGREE = abs_path('pedigree.txt')
    __PEDIGREE_POSITIONS = abs_path('pedigree.pos')
    
    '''pedfiddler's base example with an additional node with a single parent.'''
    __PEDIGREE_SINGLE_CHILD = abs_path('pedigree_single_child.txt')
    
    '''Unbalanced tree example.'''
    __PEDIGREE_UNBALANCED = abs_path('pedigree_unbalanced.txt')

    '''Unbalanced tree example.'''
    __PEDIGREE_1298 = abs_path('pedigree_1298.txt')
    __PEDIGREE_1298_GEN3 = abs_path('pedigree_1298_gen3.txt')
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_init_pedigree(self):
        '''Test reading a pedigree from file.'''
        ped = fp.read_pedigree(TestFamPlot.__PEDIGREE)
        assert_equal(ped.graph.number_of_nodes(), 22, 'Wrong number of nodes')
        assert_equal(ped.graph.number_of_edges(), 24, 'Wrong number of edges')
        assert_equal(ped.person_graph.number_of_nodes(), 15, 'Wrong number of nodes')
        assert_equal(ped.person_graph.number_of_edges(), 20, 'Wrong number of edges')

        ped = fp.read_pedigree(TestFamPlot.__PEDIGREE_SINGLE_CHILD)
        assert_equal(ped.graph.number_of_nodes(), 25, 'Wrong number of nodes')
        assert_equal(ped.graph.number_of_edges(), 27, 'Wrong number of edges')
        assert_equal(ped.person_graph.number_of_nodes(), 17, 'Wrong number of nodes')
        assert_equal(ped.person_graph.number_of_edges(), 22, 'Wrong number of edges')

    def test_generation_number(self):
        '''Test calculating generation numbers..'''
        ped = fp.read_pedigree(TestFamPlot.__PEDIGREE_UNBALANCED)
        ped_info = fp.PedigreeInfo(ped)
        cp = fp.CoordParams()
        cp.algorithm = 'default'
        fp.compute_coords(ped_info, cp)
        assert_equal(ped_info.gen_dict, {1: 0, 2: 0, 3: 1, 4: 1, 5: 1, 6: 2, 7: 2, 8: 3}, 'Wrong generation numbers')

    def test_draw_pedigree_hardcoded(self):
        '''Test drawing a pedigree as an EPS file using hard-coded positions.'''
        coord_params = fp.CoordParams()
        coord_params.algorithm = 'hardcoded'
        coord_params.coord_file = TestFamPlot.__PEDIGREE_POSITIONS
        coord_params.balance_marriages = False
        
        out = tempfile.NamedTemporaryFile(suffix='.eps', delete=False)
        ped_info = fp.draw_pedigree(TestFamPlot.__PEDIGREE, coord_params,
                                    fp.PlotParams(), fp.DEFAULT_COLORS, out)
        
        expected_coord = TestFamPlot.__read_coord_dict(coord_params.coord_file)
        test_util.assert_positions_almost_equal(ped_info.coord_dict, expected_coord,
                                                decimal=10, err_msg='Wrong coordinates')
        
    def test_draw_pedigree_default(self):
        '''Test drawing a pedigree as an EPS file using the default node placement.'''
        self.__test_draw_pedigree_default(TestFamPlot.__PEDIGREE)

    def test_draw_pedigree_1298(self):
        '''Test drawing a pedigree as an EPS file using the default node placement.'''
        self.__test_draw_pedigree_default(TestFamPlot.__PEDIGREE_1298)

    def test_draw_pedigree_1298_gen3(self):
        '''Test drawing a pedigree as an EPS file using the default node placement.'''
        self.__test_draw_pedigree_default(TestFamPlot.__PEDIGREE_1298_GEN3)

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def __read_coord_dict(in_file):
        '''Read coordinates from text file.'''
        return dict((int(row[0]), (row[1], row[2])) for row in np.genfromtxt(in_file))

    def __test_draw_pedigree_default(self, infile):
        '''Test drawing a pedigree as an EPS file using the default node placement.'''
        coord_params = fp.CoordParams()
        coord_params.algorithm = 'default'

        out = tempfile.NamedTemporaryFile(suffix='.eps', delete=True)
        ped_info = fp.draw_pedigree(infile, coord_params,
                                fp.PlotParams(), fp.DEFAULT_COLORS, out)
        #expected_coord = TestFamPlot.__read_coord_dict(TestFamPlot.__PEDIGREE_POSITIONS)
        # Positions are not unique; they depend on ordering person lists within the default algorithm
#        test_util.assert_positions_almost_equal(ped_info.coord_dict, expected_coord, 
#                                                decimal=10, err_msg='Wrong coordinates')
#        print '==============================================='
#        print ped_info.pprint()
#        print '==============================================='
        for info in ped_info.data.itervalues():
            self.assertTrue(info.coord[0] >= 0 and info.coord[1] >= 0, 
                            'Unset coordinate at node %s' % (info.name,))
        return ped_info 
    