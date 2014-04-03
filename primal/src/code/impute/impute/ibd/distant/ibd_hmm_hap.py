'''
============================================================
IBD segment identification between haplotypes using
a Hidden Markov Model (HMM).

TODO: reuse Germline to find IBD segments faster for a group
of individuals.

Created on January 25, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, hmmt, itemutil, util, matplotlib.pylab as P
from utilities.math.index_util import first_occurrence_index_byte
from impute.data.constants import MISSING

# from matplotlib import rc
# rc('text', usetex=True)

#---------------------------------------------
# Methods
#---------------------------------------------       
def prob_ibd_hmm_from_raw_input(h):
    '''A utility method that conveniently exposes the ProbIbdHmmHapCalculator class but relies on
    a minimal, serializable IbdProblem object ''h'' that encapsulates the raw arrays required
    to calculate IBD segments for a haplotype pair.'''
    return ProbIbdHmmHapCalculator(h.lam, h.f, h.x, h.p, h.h1, h.h2, snps=h.snps, e=h.params.error_rate, debug=h.params.debug).prob()

def prob_ibd_hmm(problem, hap1, hap2, snps, params):
    '''Is the segment s (defined by the SNP array snps, which may be a frame within the actual segment)
    an IBD segment between haplotypes hap1=(id1,a1) and hap2=(id2,a2) or not? Outputs an IBD probability
    estimate. The estimate is based on an HMM. Assuming the SNPs in the ''snps' array are statistically
    independent.'''
    return prob_ibd_hmm_calculator(IbdProblem(problem, hap1, hap2, snps, params)).prob()

def prob_ibd_hmm_calculator(problem, hap1, hap2, snps, params):
    '''A utility method that conveniently exposes the ProbIbdHmmHapCalculator class for debugging and plotting
    IBD posteriors in client code.'''
    h = IbdProblem(problem, hap1, hap2, snps, params)
    return ProbIbdHmmHapCalculator(h.lam, h.f, h.x, h.p, h.h1, h.h2, snps=h.snps, e=h.params.error_rate, debug=h.params.debug)

'''Kinship coefficient, given identity coefficients.'''
kinship = lambda d: d[0] + 0.5 * (d[2] + d[4] + d[6]) + 0.25 * d[7]
'''Row-wise kinship for an array of deltas.'''
kinship_vector = lambda d: d[:, 0] + 0.5 * (d[:, 2] + d[:, 4] + d[:, 6]) + 0.25 * d[:, 7]

'''Probability of IBD >= 1, given identity coefficients.'''
p_ibd = lambda d: d[0] + d[2] + d[4] + d[6] + d[7]

#---------------------------------------------
# Private Methods
#---------------------------------------------
def _error_model(e):
    '''P(O1,G2|G1,G2), where O1,O2=observed genotype, G1,G2=true genotypes, and e=genotype error rate.'''
    u = 1 - e
    e2 = e * e
    u2 = u * u
    ue = u * e
    return np.array([[u2, ue, ue, e2],
                     [ue, u2, e2, ue],
                     [ue, e2, u2, ue],
                     [e2, ue, ue, u2]])

def _emission_probability(p):
    '''P(G|S) where G=genotype state and S=identity state. Depends on the 1-allele frequency p.'''
    q = 1 - p
    pq = p * q
    return np.array([[p * p, pq, pq, q * q],
                     [p, 0, 0, q]])

def _largest_frame(problem, data_exist, debug=False):
    '''Find Largest intersection of an independent SNP frame and the data_exist array.'''
    # Fetch frames of the chromosome in question. Assuming all SNPs are on the same chromosome.
    where_data_exist = np.where(data_exist)[0]
    frames = problem.frames[problem.info.snp['chrom'][0]]
    pair_frames = [np.intersect1d(frame, where_data_exist) for frame in frames]
    frame_number = np.argmax(np.array([len(x) for x in pair_frames]))
    frame_size = len(frame)
    # Add the first and last SNP in the domain that have i,j data for full frame coverage
    start = first_occurrence_index_byte(data_exist, True, 0, 1)
    stop = first_occurrence_index_byte(data_exist, True, len(data_exist) - 1, -1) 
    frame = np.sort(np.union1d([start, stop], pair_frames[frame_number]))
    if debug:
        print 'Frame lengths', [len(x) for x in frames]
        print 'Frame #', frame_number, ', size', frame_size, 'start', start, 'stop', stop, frame.tolist()
    return frame

####################################################################################
class IbdProblem(object):
    '''A minimal serializable problem object that encapsulates the raw arrays required
    to calculate IBD segments for a haplotype pair.'''
    def __init__(self, problem, hap1, hap2, snps, params):
        '''Initialize an IBD hap problem object for the hap pair hap1=(id1, a1), hap2=(id2, a2)
        using the largest set of independent SNPs at which both haps are called, and the
        parameter struct ''params''.'''
        if params.debug: print 'Creating IbdProblem', hap1, hap2
        self.hap1, self.hap2, (id1, a1), (id2, a2), sample_id, self.params = hap1, hap2, hap1, hap2, problem.pedigree.sample_id, params
        all_h1, all_h2 = problem.h[:, id1, a1], problem.h[:, id2, a2]
        self.snps, self.num_snps = snps if snps is not None else _largest_frame(problem, (all_h1 != MISSING) & (all_h2 != MISSING), debug=params.debug), problem.num_snps
        self.h1, self.h2 = problem.h[self.snps, id1, a1], problem.h[self.snps, id2, a2]
        self.f = params.kinship(sample_id[id1], sample_id[id2])
        self.lam = problem.lam_eval(self.f)  # Recombination rate
        self.x = problem.info.snp['dist_cm'][self.snps]  # Genetic location of markers [cM] of SNPs
        self.p = problem.info.allele_frequency(1)[self.snps]  # Minor allele frequencies of SNPs
        self.bp, self.cm = problem.info.snp['base_pair'], problem.info.snp['dist_cm']
        self.d = (all_h1 == all_h2)
        # For debugging only
        if params.debug:
            self.all_h1, self.all_h2 = all_h1, all_h2

    def __repr__(self):
        return 'IbdProblem[(%d,%d), (%d,%d)]' % (self.hap1[0], self.hap1[1], self.hap2[0], self.hap2[1])

####################################################################################
class ProbIbdHmmHapCalculator(object):
    '''Calculates posterior IBD probabilities P(IBD|G(id1),G(id2) on a frame of snps)
    at each frame snp.'''
    #---------------------------------------------
    # Constants
    #---------------------------------------------

    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, lam, f, x, p, h1, h2, e, Pi=None, debug=False, snps=None):
        self.f = f
        self.x = x
        self.lam = lam
        F = np.array([1 - f, f])
        self.D = np.dot(np.ones(2)[:, np.newaxis], F[:, np.newaxis].transpose())
        self.I_minus_D = np.eye(2) - self.D
        self.exp_lam_x = np.exp(-lam * np.diff(x))
        self.p = p
        self.e = e
        self.E = _error_model(e)
        self.h1 = h1
        self.h2 = h2
        self.Pi = Pi if Pi is not None else F
        self.debug = debug
        self.snps = snps if snps is not None else np.arange(len(x))
        
        # Define HMM
        self.m = hmmt.HMM(2, A=lambda k: self.__transition_probability(k),
                          B=lambda k: self.__emission_probability(k), Pi=self.Pi, V=np.arange(4))
        self.Obs = ProbIbdHmmHapCalculator._hap_to_obs(self.h1, self.h2)
        # Output fields
        self.Gamma = None
        self.p_ibd_gamma = None
        self.Q_star = None
        self.p_ibd = None
        
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __repr__(self): return 'ProbIbdHmmHapCalculator[lam=%.2f, Delta=%s]' % (self.lam, repr(self.Delta))

    @property
    def data(self): return np.array([self.x, np.concatenate((self.exp_lam_x, [0])), self.p]).transpose()

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def print_table(self):
        '''Print a table of probabilities at each SNP.'''
        options = np.get_printoptions()
        np.set_printoptions(precision=3, suppress=True, threshold=np.nan, linewidth=200)
        print 'lambda = %.2f, f = %.2f, eps = %.1e' % (self.lam, self.f, self.e)
        print 'Viterbi path (frame SNPs): ' + ' -> '.join(map(lambda x: '%d (%d-%d)' % (x[0], x[1][0], x[1][1]),
                                                  itemutil.groupby_with_range(self.Q_star)))
        print 'Viterbi path (SNPs):       ' + ' -> '.join(map(lambda x: '%d (%d-%d)' % (x[0], self.snps[x[1][0]], self.snps[x[1][1]]),
                                                  itemutil.groupby_with_range(self.Q_star)))
        print '      %-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s' % \
        ('t', 'SNP#', 'H1', 'H2', 'x', 'lam*dx', 'p', 'Gam0', 'Gam1', 'p(IBD)', 'Viterbi', 'IBD?')
        print np.concatenate((np.arange(len(self.x))[np.newaxis].transpose(),
                              self.snps[np.newaxis].transpose(),
                              self.Obs[np.newaxis].transpose() / 2,
                              self.Obs[np.newaxis].transpose() % 2,
                              self.x[np.newaxis].transpose(),
                              np.concatenate((self.lam * np.diff(self.x), [0]))[np.newaxis].transpose(),
                              self.p[np.newaxis].transpose(),
                              self.Gamma.transpose(),
                              self.p_ibd_gamma[np.newaxis].transpose(),
                              (self.Q_star)[np.newaxis].transpose(),
                              self.p_ibd_viterbi[np.newaxis].transpose()
                              ), axis=1)
        util.set_printoptions(options)
        
    def plot(self):
        '''Plot the IBD probabilities at each SNP.'''
        P.clf()
        X = self.snps
        for i in xrange(2):
            P.plot(X, self.Gamma[i, :], label=r'$\gamma_%d$' % (i + 1,))
        P.plot(X, self.p_ibd_gamma, 'k-', linewidth=2, label='P(IBD)')
        P.plot(X, (self.Q_star + 1) / 9., 'b-', linewidth=2, label='Viterbi')
        P.ylim([-0.1, 1.1])
        P.xlabel('Genetic Distance (x)')
        P.ylabel('State Probability')
        P.title('HMM Posterior Probabilities')
        P.legend(loc='lower right')

    def prob(self):
        '''Return the posterior IBD probability at each marker. First, we find the IBD SNPs
        using the Viterbi path; at each of them, we return the appropriate Gamma sum as the 
        posterior; outside segments, the posterior is set to 0. 
        
        This confounds two HMM optimality definitions. The reason for doing that is that Viterbi seems
        more robust to finding IBD -segments- (no need to threshold it and obtain fragmented segments
        like the Gamma posterior would yield); but Viterbi gives no confidence intervals, so we use
        Gamma once we are in a Viterbi-flagged segment.'''
        p_ibd = self.__calculate_prob_gamma()
        # Zero out the probability of IBD outside segments flagged by Viterbi
        is_ibd_viterbi = self.__calculate_prob_viterbi()
        # Print debugging table of HMM information
        if self.debug:
            self.print_table()
            # P.figure(100)
            # self.plot()
        p_ibd[~is_ibd_viterbi] = 0.0
        return p_ibd
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    '''Convert from 1/2-recoded haplotypes to observation codes.''' 
    @staticmethod
    def _hap_to_obs(h1, h2): return 2 * h1 + h2 - 3
    def __transition_probability(self, k): return self.D + self.I_minus_D * self.exp_lam_x[k]
    def __emission_probability(self, k): return np.dot(_emission_probability(self.p[k]), self.E)
    
    def __calculate_prob_viterbi(self):
        '''Return the posterior IBD probability at each marker using the Viterbi path.
        Here it is 1 if the state is an IBD>=1 state, otherwise 0. Sets the self.Q_star field.'''
        # Remove observations that are IBS=0 - debugging
        Q_star, _, _ = hmmt.viterbi(self.m, self.Obs, scaling=1)
        self.Q_star = Q_star
        self.p_ibd_viterbi = (Q_star == 1)
        return self.p_ibd_viterbi
     
    def __calculate_prob_gamma(self):
        '''Return the posterior IBD probability at each marker using individual-state
        posterior probabilities (Gamma). Sets the self.Gamma field.'''
        # Remove observations that are IBS=0 - debugging
        _, Alpha, c = hmmt.forward(self.m, self.Obs, scaling=1)
        Beta = hmmt.backward(self.m, self.Obs, c)
        Gamma, _ = hmmt.individually_optimal_states(Alpha, Beta)
        self.Gamma = Gamma
        self.p_ibd_gamma = Gamma[1, :]
        return self.p_ibd_gamma 
