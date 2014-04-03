'''
============================================================
Algorithm parameters. Used by all high-level packages
(ibd, phasing, etc.).

Created on September 14, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, util, db_gene

class PhaseParam(util.Struct):
    '''A struct of parameters of the phasing algorithm.'''
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    class TEMPLATE_MODE: 
        '''Sib comparison mode: ALL=use all template children (expensive but more accurate);
        ONE=only the most-filled one (faster but risky); TWO=the two most-filled (a compromise).'''
        ONE, TWO, ALL = range(3)

    '''Default fill threshold of parent [het snp] haplotypes.'''
    HET_FILL_THRESHOLD = 0.5 
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, **entries):
        '''Set the default parameters.'''
        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        # Debugging
        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        '''Run a single stage or all stages (if 0)'''
        self.stage = 0
        
        '''Stage Timing printout flags'''
        self.print_times = False
        
        '''Debugging printouts flag'''
        self.debug = False
        
        '''Process only this family.'''
        self.single_family = None
        
        '''Process only this member or its families.'''
        self.single_member = None
        
        '''If not None, phase only the samples in this array, assuming the rest are already
        phased, and their haplotypes are stored in the input genotype file.'''
        self.selected_samples = None
        
        '''An isolated SNP index to debug.'''
        # self.snp_test_index = -1
        
        '''Number of parallel processes to spawn.'''
        self.num_processes = 1

        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        # IBD Segments and Fill
        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        '''In sib comparison: use all template children (expensive but more accurate),
        only the most-filled one (faster but risky), or the two most-filled (a compromise)'''
        self.template_mode = PhaseParam.TEMPLATE_MODE.ALL

        '''Expected genotype error rate.'''
        self.error_rate = 0.01
        
        '''Minimum haplotype fraction at het SNPs for a segment to be called phased.'''
        self.het_fill_threshold = 0.5

        '''Target fill at samples to be phased using distant relatives.'''
        self.target_fill = 0.95
        
        '''Minimum haplotype fraction (at all SNPs) required to consider a surrogate parent for IBD calculations.'''
        self.surrogate_parent_fill_threshold = 0.9
        
        '''A list of (surrogate_parent_fill_threshold, max_path_length, target_fill) tuples,
        each corresponding to distant phasing pass.'''
        self.distant_phasing_params = [(0.9, 6, 0.95), (0.92, 7, 0.92)]
                        
        '''Minimum # meioses to search for surrogate parents whose IBD segments cover the sample to be phased.''' 
        self.min_path_length = 3

        '''Maximum # meioses to search for surrogate parents  whose IBD segments cover the sample to be phased.''' 
        self.max_path_length = 6
        
        '''Minimum probability of IBD to call an IBD segment.'''
        self.prob_ibd_threshold = 0.95
        
        '''Minimum IBS segment size to consider as IBD segment.'''
        self.min_ibs_len_snp = 400

        '''IBD segment length unit'''
        self.len_unit = 'mbp'  # 'cm'

        '''Minimum IBD segment length to look for between distant relatives [centi-Morgans]'''
        self.min_len = 1.0

        '''Minimum IBD sub-segment length to look for within a segment identified by HMM [centi-Morgans].
        As of 5-MAR-13, 1.0/0.4 segment thresholds are optimal, at least for chr22.'''
        self.min_len_sub_segment = 0.4 * self.min_len
        
        '''Minimum segment length to look for between siblings [mega base pairs]'''
        self.min_segment_length_sibs = 5.0
        
        '''Minimum number of haplotypes required to agree on an IBD concensus haplotype.'''
        self.min_consensus_samples = 3

        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        # GERMLINE IBD Segments
        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        '''Initial GERMLINE slice length (I).
        * Chance of obtaining a false segment in the initial slice among N samples: 
                            1-(1-2**-I)**(N*(N-1)/2) (e.g., I=20, N=10: 4.3e-5)
        * Chance of missing a segment in the initial slice among N samples with genotype error rate e: 
                            1-(1-e)**(2*I*N) (e.g., I=20, N=10, e=0.01: .86. But at least a -short- segment is missing only.)
        '''
        self.initial_slice_size = 20
        self.slice_size = 100

        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        # PARENT-OF-ORIGIN ALIGNMENT
        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        '''POO measure used during phasing stage 4 to align quasi-founder sibs. Based on
        paternal haplotype coverage. Should be lower than the SNP measure threshold below
        due to the difference in measure definitions.'''
        self.poo_coloring_measure_threshold = 0.3

        '''POO SNP measure threshold above which we set phase=1 and below -(which) we set phase=-1.'''
        self.poo_snp_measure_threshold = 0.5
        '''POO measure threshold above which we include SNPs in POO determination. In [0,1]. 1 means
        perfect alignment of a asample, 0=no alignment of a sample.'''
        self.poo_threshold = 0.75
        '''SNP step size of down-sampling all chromosome SNPs.'''
        self.poo_snp_step_size = 1
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        # External Data
        # %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        '''Path to identity coefficient file.'''
        self.id_coef_file = os.environ['OBER_OUT'] + '/phasing/hutt.id'

        '''Path to kinship file.'''
        self.kinship_file = os.environ['OBER_OUT'] + '/phasing/hutt.kinship'
        
        util.Struct.__init__(self, **entries)
        
        # Lazily-initialized dependent properties
        self._idcoef_dao = None
        self._kinship_dao = None
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    @property
    def selected_mode(self):
        '''Is phasing selected samples mode on?'''  
        return self.selected_samples is not None

    def id_coefs(self, id1, id2):
        '''Load the condensed identity coefficients (lambda, (Delta1,...,Delta9)) between id1 and id2.'''
        if not self._idcoef_dao:
            self._idcoef_dao = db_gene.snp.file_dao.IdCoefDao(self.id_coef_file)
        return self._idcoef_dao.id_coefs(id1, id2) 

    def kinship(self, id1, id2):
        '''Load the kinship coefficients between id1 and id2.'''
        if not self._kinship_dao:
            self._kinship_dao = db_gene.snp.file_dao.KinshipDao(self.kinship_file)
        return self._kinship_dao.kinship(id1, id2) 
