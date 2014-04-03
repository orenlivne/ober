'''
============================================================
An extension of the hmm module to HMM with time-dependent
transition and emission matrices.

Created on January 5, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from hmm import symbol_index
import numpy
# from matplotlib import rc
# rc('text', usetex=True)

class HMM:
    """
    Creates and maintains a hidden Markov model.  This version assumes the every state can be
    reached DIRECTLY from any other state (ergodic).  This, of course, excludes the start state.
    Hence the state transition matrix, A, must be N X N .  The observable symbol probability
    distributions are represented by an N X M matrix where M is the number of observation
    symbols.  

                  |a_11 a_12 ... a_1N|                   |b_11 b_12 ... b_1M|   
                  |a_21 a_22 ... a_2N|                   |b_21 b_22 ... b_2M|
              A = | .    .        .  |               B = | .    .        .  |
                  | .         .   .  |                   | .         .   .  |
                  |a_N1 a_N2 ... a_NN|                   |b_N1 b_N2 ... b_NM|
          
           a_ij = P(q_t = S_j|q_t-1 = S_i)       b_ik = P(v_k at t|q_t = S_i)
        where q_t is state at time t and v_k is k_th symbol of observation sequence

    Here, A and B are also functions of the state (A(t): transition from time t to t+1; B(t):
    emission at time t.  
    """
    def __init__(self, n_states=1, **args):
        """
        :Keywords:
          - `n_states` - number of hidden states
          - `V` - list of all observable symbols
          - `A` - transition matrix
          - `B` - observable symbol probability distribution
          - `D` - dimensionality of continuous observations
          """
    
        self.N = n_states  # Number of hidden states

        # Initialize observable symbol set parameters
        self.V = args[ 'V' ] 
        self.M = len(self.V)
        self.symbol_map = dict(zip (self.V, range(len(self.V))))
            
            
        # Initialize transition probability matrix
        if 'A' in args:
            self.A = args[ 'A' ] 
            # Test that A is a function and that A(t=0) (for instance) has the right shape
            # assert hasattr(self.A, '__call__') 
            assert numpy.shape(self.A(0)) == (self.N, self.N)
        else:
            # Randomly initialize matrix and normalize so sum over a row = 1. Constant in time.
            raw_A = numpy.random.uniform(size=self.N * self.N).reshape((self.N, self.N))
            const_A = (raw_A.T / raw_A.T.sum(0)).T  
            if n_states == 1:
                const_A.reshape((1, 1))
            self.A = lambda t: const_A
            
        # Initialize observable symbol probability distributions
        if 'B' in args:
            self.B = args[ 'B' ]
            if n_states > 1:
                # assert hasattr(self.B, '__call__') 
                assert numpy.shape(self.B(0)) == (self.N, self.M)
            else:
                assert numpy.shape(self.B(0)) == (1, self.M)
                # self.B = numpy.reshape(self.B, (1, self.M))
        else:
            # initialize distribution. Constant in time.
            B_raw = numpy.random.uniform(0, 1, self.N * self.M).reshape((self.N, self.M))
            const_B = (B_raw.T / B_raw.T.sum(0)).T
            self.B = lambda t: const_B

        # Initialize the intitial state distribution
        if 'Pi' in args:
            self.Pi = args[ 'Pi' ]
            assert len(self.Pi) == self.N
        else:
            # initialize to uniform distribution
            self.Pi = numpy.array (1.0 / self.N).repeat(self.N)

        if 'Labels' in args:
            self.Labels = args[ 'Labels' ]
        else:
            self.Labels = range(self.N)

    def __repr__(self):
        retn = ""
        retn += "num hiddens: %d\n" % (self.N) + \
                "symbols: %s\n" % (self.V) + \
                "\nA(t=0):\n %s\n" % (str(self.A(0))) + \
                "Pi:\n %s" % (str(self.Pi))
        
        return retn

def forward(hmm, Obs, scaling=True):
    """
    Calculate the probability of an observation sequence, Obs,
    given the model, P(Obs|hmm).
    Obs: observation sequence
    hmm: model
    returns: P(Obs|hmm)
    """
    # Number of states in observation sequence
    T = len(Obs)
    # Get index sequence of observation sequence to access the observable symbol prob distribution matrix
    Obs = symbol_index(hmm, Obs)
    # Allocate output array, scaling array 
    Alpha = numpy.zeros([ hmm.N, T ], float)
    if scaling:
        c = numpy.zeros([ T ], float)
   
#    print 'forward sweep'
#    print 'Obs', Obs
    for t in xrange(T):
#        print 't', t
#        if t > 0:
#            print 'A', hmm.A(t - 1)
#            print 'B', hmm.B(t)
#            print 'Obs[t]', Obs[t]
#            print 'B[:, Obs[t]]', hmm.B(t)[ :, Obs[ t ] ]
        Alpha[ :, t ] = (hmm.Pi if t == 0  # Base Case
                         else numpy.dot(Alpha[ :, t - 1 ], hmm.A(t - 1))  # Induction Step:
                         ) * hmm.B(t)[ :, Obs[ t ] ]
        if scaling:
            c[ t ] = 1.0 / numpy.sum(Alpha[ :, t ])
            Alpha[ :, t] = Alpha[ :, t] * c[ t ]
            #print t, Alpha[:, t], c[t]
            
    return (-(numpy.sum(numpy.log(c))), Alpha, c) if scaling else (numpy.sum(Alpha[ :, T - 1 ]), Alpha)
            

def backward(hmm, Obs, c=None):
    """
    Calculate the probability of a partial observation sequence
    from t+1 to T, given some state t.
    Obs: observation sequence
    hmm: model
    c: the scaling coefficients from forward algorithm
    returns: B_t(i) 
    """
    # Number of states in observation sequence
    T = len(Obs)
    # Get index sequence of observation sequence to access the observable symbol prob distribution matrix
    Obs = symbol_index(hmm, Obs)
    # Allocate output array
    Beta = numpy.zeros([ hmm.N, T ], float) 

    # Inductive Step:
    for t in reversed(xrange(T)):
        # Base Case
        Beta[ :, t ] = 1.0 if t == T - 1 else \
        numpy.dot(hmm.A(t), (hmm.B(t + 1)[ :, Obs[ t + 1 ] ] * Beta[ :, t + 1 ]))  # Inductive Step
        if c is not None:
            # Rarely, Beta might blow up when Alpha is exactly zero. This is not a problem without
            # scaling. For simplicity, cap Beta.
            # TODO: it might be safer to restart the algorithm with no scaling instead.  
            Beta[ :, t ] = numpy.minimum(Beta[ :, t ] * c[ t ], 1e+100) 
#            print t, Beta[:, t], c[t]
            #pass
            
    return Beta
            
def individually_optimal_states(Alpha, Beta):
    '''Compute the most probable states given an observation sequence, and the
    state probability matrix Gamma.'''        
    Gamma_raw = Alpha * Beta
    return (Gamma_raw / Gamma_raw.sum(0), numpy.argmax(Gamma_raw, 0))

def viterbi(hmm, Obs, scaling=True):
    """
    Calculate P(Q|Obs, hmm) and yield the state sequence Q* that
    maximizes this probability. 
    Obs: observation sequence
    hmm: model
    """
    # Number of states in observation sequence
    T = len(Obs)
    # Get index sequence of observation sequence to access the observable symbol prob distribution matrix
    Obs = symbol_index(hmm, Obs)
    
    # Initialization
    # Delta[ i,j ] = max_q1,q2,...,qt P( q1, q2,...,qt = i, O_1, O_2,...,O_t|hmm )
    # this is the highest prob along a single path at time t ending in state S_i
    Delta = numpy.zeros([ hmm.N, T ], float)
    # Track Maximal States
    Psi = numpy.zeros([ hmm.N, T ], int)

#    print 'scaling', scaling
    if scaling:
        Delta[ :, 0 ] = _log(hmm.Pi) + _log(hmm.B(0)[ :, Obs[ 0 ] ])
    else:
        Delta[ :, 0 ] = hmm.Pi * hmm.B(0)[ :, Obs[ 0] ]
#    print 't', 0
#    print 10 ** Delta[:, 0][numpy.newaxis].transpose()
    
    # Inductive Step:
    if scaling:
        for t in xrange(1, T):
            nus = Delta[ :, t - 1 ] + _log(hmm.A(t - 1).transpose())
            Delta[ :, t ] = nus.max(1) + _log(hmm.B(t)[ :, Obs[ t ] ])
            Psi[ :, t ] = nus.argmax(1)
#            print 't', t
#            print 'A',
#            print hmm.A(t - 1)
#            print 'nus'
#            print 10 ** nus
#            print 'Delta'
#            print 10 ** Delta[:, t][numpy.newaxis].transpose()
#            print 'Psi'
#            print Psi[:, t][numpy.newaxis].transpose()
    else:
        for t in xrange(1, T):
            nus = Delta[ :, t - 1 ] * hmm.A(t - 1).transpose()
            Delta[ :, t ] = nus.max(1) * hmm.B(t)[ :, Obs[ t ] ]
            Psi[ :, t ] = nus.argmax(1)
        
#    print 10 ** Delta
#    print Psi
    
    # Calculate State Sequence, Q*:
    Q_star = [ numpy.argmax(Delta[ :, T - 1 ]) ]
    for t in reversed(xrange(T - 1)) :
        Q_star.insert(0, Psi[ Q_star[ 0 ], t + 1 ])

    return (numpy.array(Q_star), Delta, Psi)

def _log(A):
    return numpy.log10(A + 1e-15)

def log_likelihood(hmm, Obs, q):
    '''Return log P(Obs|q) (if interpreted as a function of q for fixed Obs, this is the likelihood of the state sequence q.'''
    p = _log(hmm.Pi[q[0]]) + _log(hmm.B(0)[q[0], Obs[0]])
    for t in xrange(1, len(Obs)):
        p += _log(hmm.A(t - 1)[q[t - 1], q[t]])
        p += _log(hmm.B(t))[q[t], Obs[t]]
    return p
    
