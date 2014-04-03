'''
============================================================
Impute project - main test suite.

Created on Jun 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import test_util as tu
from tests.data.pedigree.TestPedigree import TestPedigree
from tests.data.pedigree.TestPedigreeTools import TestPedigreeTools
from tests.data.pedigree.TestPedigreePlotLaplacian import TestPedigreePlotLaplacian
from tests.data.genotype.TestGenotypeTools import TestGenotypeTools
from tests.data.genotype.TestGenotype import TestGenotype
from tests.data.problem.TestProblem import TestProblem
from tests.data.TestIo import TestIo
from tests.ibd.TestSegment import TestSegment
from tests.ibd.TestIbd import TestIbd
from tests.ibd.TestIbdIbdld import TestIbdIbdld
from tests.ibd.TestIbdGermline import TestIbdGermline
from tests.ibd.TestIbdGermline1 import TestIbdGermline1
from tests.ibd.TestIbdGermline1Restricted import TestIbdGermline1Restricted
from tests.ibd.TestIbdGermline2Small import TestIbdGermline2Small
from tests.ibd.distant.TestIbdDistantFamily import TestIbdDistantFamily
from tests.ibd.distant.TestIbdHmm import TestIbdHmm
from tests.ibd.distant.TestIbdHmmHap import TestIbdHmmHap
from tests.ibd.distant.TestIbdDistantHapFamily import TestIbdDistantHapFamily
from tests.phasing.TestPhasePrepare import TestPhasePrepare
from tests.phasing.trivial.TestPhaseTrivialTrio import TestPhaseTrivialTrio
from tests.phasing.trivial.TestPhaseTrivialDuo import TestPhaseTrivialDuo
from tests.phasing.trivial.TestPhaseTrivialAll import TestPhaseTrivialAll
from tests.phasing.family.TestPhaseFamilyParentChild import TestPhaseFamilyParentChild
from tests.phasing.family.TestPhaseFamily import TestPhaseFamily
from tests.phasing.family.TestPhaseFamilyCompareChild import TestPhaseFamilyCompareChild
from tests.phasing.family.TestPhaseFamilyCompareChild0 import TestPhaseFamilyCompareChild0
from tests.phasing.family.TestPhaseFamilyCompareSibs import TestPhaseFamilyCompareSibs
from tests.phasing.distant.TestPhaseDistant import TestPhaseDistant
from tests.phasing.distant.TestPhaseDistantFounderSibs import TestPhaseDistantFounderSibs
from tests.phasing.TestPhasePipeline import TestPhasePipeline
from tests.validation.TestValidation import TestValidation
from tests.ibd.TestSmartSegmentSet import TestSmartSegmentSet
from tests.ibd.index.TestSegmentIndex import TestSegmentIndex
from tests.imputation.TestImputeStats import TestImputeStats
from tests.imputation.TestImputation import TestImputation
from tests.ibd.index.TestSmartSegmentSetCython import TestSmartSegmentSetCython
from tests.color.TestHapColorGrouping import TestHapColorGrouping
from tests.color.TestHapColor import TestHapColoring
from tests.phasing.distant.TestPhaseDistantFounderSibs1049 import TestPhaseDistantFounderSibs1049
from tests.cgi.TestOverrideImputedByPlink import TestOverrideImputedByPlink
# from tests.imputation.TestImputationSnp import TestImputationSnp

def suite():
    return tu.load_tests_from_classes([TestPedigree, TestPedigreeTools,
                                       TestPedigreePlotLaplacian,
                                       TestGenotype, TestGenotypeTools,
                                       TestIo, TestProblem,
                                       
                                       TestSegment, TestHapColoring, TestHapColorGrouping,
                                       TestIbd, TestIbdDistantFamily, TestIbdDistantHapFamily,
                                       TestIbdHmm, TestIbdHmmHap, TestIbdIbdld,
                                       TestIbdGermline, TestIbdGermline1, TestIbdGermline1Restricted, TestIbdGermline2Small,

                                       TestPhasePrepare,
                                       TestPhaseTrivialDuo, TestPhaseTrivialTrio,
                                       TestPhaseTrivialAll,
                                       TestPhaseFamilyParentChild,
                                       TestPhaseFamily, TestPhaseFamilyCompareChild, TestPhaseFamilyCompareChild0,
                                       TestPhaseFamilyCompareSibs,
                                       TestPhaseDistant, TestPhaseDistantFounderSibs, TestPhaseDistantFounderSibs1049,
                                       TestPhasePipeline,
                                       
                                       TestSmartSegmentSet,
                                       TestSegmentIndex, TestSmartSegmentSetCython,
                                       
                                       TestValidation,
                                       TestImputation,
#                                        TestImputationSnp, 
                                       TestImputeStats,
                                       #
                                       TestOverrideImputedByPlink
                                       ])
