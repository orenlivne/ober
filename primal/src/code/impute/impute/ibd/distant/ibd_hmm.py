'''
============================================================
IBD segment identification between genotypes using
a Hidden Markov Model (HMM).

TODO: reuse Germline to find IBD segments faster for a group
of individuals.

Created on September 14, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, hmmt, itemutil, util, matplotlib.pylab as P, itertools
from impute.tools import recode
from collections import OrderedDict

#---------------------------------------------
# Methods
#---------------------------------------------
####################################################################################
def prob_ibd_hmm(problem, id1, id2, snps, params):
    '''Is the segment s (defined by the SNP array snps, which may be a frame within the actual segment)
    an IBD segment between samples id1 and id2 or not? Outputs an IBD probability estimate.
    The estimate is based on an HMM. Assuming loci are statistically independent.'''
    sample_id = problem.pedigree.sample_id
    lam, Delta = params.id_coefs(sample_id[id1], sample_id[id2])
    x = problem.info.snp['dist_cm'][snps]  # [cM]
    p = problem.info.allele_frequency(1)[snps]
    g1, g2 = problem.g[snps, id1, :], problem.g[snps, id2, :]
    calculator = ProbIbdHmmCalculator(0.4 * lam, Delta, x, p, g1, g2,
                                      e=params.error_rate, debug=params.debug, snps=snps)
    return calculator.prob()

#---------------------------------------------
# Private Methods
#---------------------------------------------
class ProbIbdHmmCalculator(object):
    '''Calculates posterior IBD probabilities P(IBD|G(id1),G(id2) on a frame of snps)
    at each frame snp.'''
    #---------------------------------------------
    # Constants
    #---------------------------------------------

    # The T-state (IBS state) corresponding to hash code of (recoded genotype A, recoded genotype B).
    # Hash code = 3*(rA-2) + rB-2 = 3*rA + rB - 8 where rA, rB are the recoded genotypes.
    # States T=2,6 correspond to IBS=0, the rest to IBS>=1.
    __HASH_TO_T_STATE = np.array([1, 3, 2, 5, 7, -5, -2, -3, -1])

    # List of HMM observation symbols corresponding to the columns of the emission probability matrix B
    __T_STATE = np.array([1, -1, 2, -2, 3, -3, 5, -5, 7])
    # Corresponding genotypes, for pretty-printouts
    _T_STATE_G = OrderedDict([(+1, (11, 11)),
                              (-1, (22, 22)),
                              (+2, (11, 22)),
                              (-2, (22, 11)),
                              (+3, (11, 12)),
                              (-3, (22, 12)),
                              (+5, (12, 11)),
                              (-5, (12, 22)),
                              (+7, (12, 12))])

    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, lam, Delta, x, p, g1, g2, e, Pi=None, debug=False, snps=None):
        self.Delta = Delta 
        self.x = x
        self.lam = lam
        self.D = np.dot(np.ones(9)[:, np.newaxis], Delta[:, np.newaxis].transpose())
        self.I_minus_D = np.eye(9) - self.D
        self.lam_x = lam * np.diff(x)
        self.p = p
        self.e = e
        self.E = emission_error(e)
        self.r1 = recode.recode_single_genotype(g1)
        self.r2 = recode.recode_single_genotype(g2)
        self.Pi = Pi if Pi is not None else Delta
        self.debug = debug
        self.snps = snps if snps is not None else np.arange(len(x))
        
        # Define HMM
        self.m = hmmt.HMM(9, A=lambda k: self.__transition_probability(k),
                          B=lambda k: self.__emission_probability(k),
                          Pi=self.Pi, V=ProbIbdHmmCalculator.__T_STATE)
        self.Obs = ProbIbdHmmCalculator.__HASH_TO_T_STATE[3 * self.r1 + self.r2 - 8]
        # Output fields
        self.Gamma = None
        self.p_ibd_gamma = None
        self.Q_star = None
        self.p_ibd = None
        
    #---------------------------------------------
    # Operators
    #---------------------------------------------
    def __repr__(self):
        s = 'ProbIbdHmmCalculator[lam=%.2f, Delta=%s]' % (self.lam, repr(self.Delta))
        return s 

    @property
    def data(self):        
        return np.array([self.x, np.concatenate((self.lam_x, [0])), self.p]).transpose()

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def remove_non_ibs(self):
        '''Remove observations that are IBS=0. For debugging only.'''
        Obs = ProbIbdHmmCalculator.__HASH_TO_T_STATE[3 * self.r1 + self.r2 - 8]
        ok = np.where((Obs != 2) & (Obs != -2))
        self.snps = self.snps[ok]
        Obs = Obs[ok]
        self.p = self.p[ok]
        self.x = self.x[ok]
        self.lam_x = self.lam * np.diff(self.x)
        self.r1 = self.r1[ok]
        self.r2 = self.r2[ok]
    
    def print_table(self):
        '''Print a table of probabilities at each SNP.'''
        options = np.get_printoptions()
        np.set_printoptions(precision=3, suppress=True, threshold=np.nan, linewidth=200)
        print 'lambda = %s, Delta = %s, eps = %.1e' % (self.lam, repr(self.Delta)[6:-1], self.e)
        print 'Viterbi path (frame SNPs): ' + ' -> '.join(map(lambda x: '%d (%d-%d)' % (x[0], x[1][0], x[1][1]),
                                                  itemutil.groupby_with_range(self.Q_star + 1)))
        print 'Viterbi path (SNPs):       ' + ' -> '.join(map(lambda x: '%d (%d-%d)' % (x[0], self.snps[x[1][0]], self.snps[x[1][1]]),
                                                  itemutil.groupby_with_range(self.Q_star + 1)))
        print '    %-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s' % \
        ('t', 'SNP#', 'Obs', 'G1', 'G2', 'lam*dx', 'p',
         'Gam1', 'Gam2', 'Gam3', 'Gam4', 'Gam5', 'Gam6', 'Gam7', 'Gam8', 'Gam9',
         'p(IBD)', 'Viterbi', 'IBD?')
        print np.concatenate((np.arange(len(self.x))[np.newaxis].transpose(),
                              self.snps[np.newaxis].transpose(),
                              self.Obs[np.newaxis].transpose(),
                              np.array([ProbIbdHmmCalculator._T_STATE_G[t][0] for t in self.Obs])[np.newaxis].transpose(),
                              np.array([ProbIbdHmmCalculator._T_STATE_G[t][1] for t in self.Obs])[np.newaxis].transpose(),
                              np.concatenate((self.lam_x, [0]))[np.newaxis].transpose(),
                              self.p[np.newaxis].transpose(),
                              self.Gamma.transpose(),
                              self.p_ibd_gamma[np.newaxis].transpose(),
                              (self.Q_star + 1)[np.newaxis].transpose(),
                              self.p_ibd_viterbi[np.newaxis].transpose()
                              ), axis=1)
        util.set_printoptions(options)
        
    def plot(self):
        '''Plot the IBD probabilities at each SNP.'''
        P.figure(2)
        P.clf()
#        from matplotlib import rc
#        rc('text', usetex=True)
        X = self.snps
        for i in xrange(9):
            if np.max(self.Gamma[i, :]) >= 0.01:
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
        is_ibd_viterbi = self.__calculate_prob_viterbi()
        p_ibd[~is_ibd_viterbi] = 0.0
        # Print debugging table of HMM information
        if self.debug:
            self.print_table()
            # P.figure(100)
            # self.plot()
        return p_ibd
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    def __transition_probability(self, k):
        return self.D + self.I_minus_D * np.exp(-self.lam_x[k])
        
    def __emission_probability(self, k):
        B = np.dot(emission_probability(self.p[k]), self.E)
        # B = emission_probability_random_error(self.p[k], self.e)
        return B
    
    def __calculate_prob_viterbi(self):
        '''Return the posterior IBD probability at each marker using the Viterbi path.
        Here it is 1 if the state is an IBD>=1 state, otherwise 0. Sets the self.Q_star field.'''
        # Remove observations that are IBS=0 - debugging
        Q_star, _, _ = hmmt.viterbi(self.m, self.Obs, scaling=1)
        self.Q_star = Q_star
        self.p_ibd_viterbi = (Q_star == 0) | (Q_star == 2) | (Q_star == 4) | (Q_star == 6) | (Q_star == 7)
        return self.p_ibd_viterbi
     
    def __calculate_prob_gamma(self):
        '''Return the posterior IBD probability at each marker using individual-state
        posterior probabilities (Gamma). Sets the self.Gamma field.'''
        # Remove observations that are IBS=0 - debugging
        _, Alpha, c = hmmt.forward(self.m, self.Obs, scaling=1)
        Beta = hmmt.backward(self.m, self.Obs, c)
        Gamma, _ = hmmt.individually_optimal_states(Alpha, Beta)
        self.Gamma = Gamma
        self.p_ibd_gamma = np.sum(Gamma[np.array([0, 2, 4, 6, 7]), :], axis=0)
        return self.p_ibd_gamma 
    
def error_model(e):
    '''P(O|G), where O=observed genotype, G=true genotype, and e=genotype error rate.'''
    u = 1 - e
    u2 = u * u
    e2 = e * e
    ue = u * e
    return np.array([[u2, 2 * ue, e2],
                     [ue, u2 + e2, ue],
                     [e2, 2 * ue, u2]])

def emission_error(e):
    '''P(O1,G2|G1,G2), where O1,O2=observed genotype, G1,G2=true genotypes, and e=genotype error rate.'''
    a = ProbIbdHmmCalculator._T_STATE_G
    k = np.array([map(lambda x: (x / 10) + (x % 10) - 2, a[x]) for x in a.iterkeys()])
    E = error_model(e)
    F = np.zeros((9, 9))
    for (i, j) in itertools.product(xrange(9), xrange(9)):
        F[i, j] = E[k[i, 0], k[j, 0]] * E[k[i, 1], k[j, 1]]
    return F

def emission_probability(p):
    '''P(G|S) where G=genotype state and S=identity state. Depends on the 1-allele frequency p.'''
    q = 1 - p
    p2 = p * p
    p3 = p2 * p
    p4 = p3 * p
    q2 = q * q
    q3 = q2 * q
    q4 = q3 * q
    pq = p * q
    return np.array([[p , q , 0      , 0       , 0         , 0          , 0         , 0         , 0          ],
                     [p2, q2, pq     , pq      , 0         , 0          , 0         , 0         , 0          ],
                     [p2, q2, 0      , 0       , pq        , pq         , 0         , 0         , 0          ],
                     [p3, q3, p * q2 , p2 * q  , 2 * p2 * q, 2 * p * q2 , 0         , 0         , 0          ],
                     [p2, q2, 0      , 0       , 0         , 0          , pq        , pq        , 0          ],
                     [p3, q3, p2 * q , p * q2  , 0         , 0          , 2 * p2 * q, 2 * p * q2, 0          ],
                     [p2, q2, 0      , 0       , 0         , 0          , 0         , 0         , 2 * pq     ],
                     [p3, q3, 0      , 0       , p2 * q    , p * q2     , p2 * q    , p * q2    , pq         ],
                     [p4, q4, p2 * q2, p2 * q2 , 2 * p3 * q, 2 * p * q3 , 2 * p3 * q, 2 * p * q3, 4 * p2 * q2]])

def emission_probability_random_error(p, e):
    q = 1 - p
    p2 = p * p
    p3 = p2 * p
    p4 = p3 * p
    q2 = q * q
    q3 = q2 * q
    q4 = q3 * q
    b9 = [p4, q4, p2 * q2, p2 * q2 , 2 * p3 * q , 2 * p * q3 , 2 * p3 * q , 2 * p * q3 , 4 * p2 * q2 ]
    random_g = np.array([b9, b9, b9, b9, b9, b9, b9, b9, b9])
    return (1 - e) * emission_probability(p) + e * random_g
