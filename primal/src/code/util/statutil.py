'''
============================================================
Statistics and numerical linear algebra and utilities.

Created on October 11, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, collections
import itertools as it
from scipy.special import betainc
from scipy import stats

#---------------------------------------------
# Frequencies and Counts
#---------------------------------------------
def group_by_value(values):
    '''Convert an array of integers to a dictionary (value,count). Perhaps not the fastest implementation; see
    http://stackoverflow.com/questions/4651683/numpy-grouping-using-itertools-groupby-perfo'''
    counts = collections.defaultdict(int)
    for v in values: counts[v] += 1
    return counts

def hist(values, n):
    '''Convert a vector 'values' of integers in the range 0..n-1 into a histogram vector of size n,
    whose ith entry is the number of times i occurs in 'values'. If values is a matrix,
    histograms are generated for each row of values.'''
    if values.ndim == 1:
        histogram = np.zeros((n,), dtype=np.uint)
        for (k, v) in group_by_value(values).iteritems():
            histogram[k] = v
        return histogram
    else:
        histogram = np.zeros((values.shape[0], n), dtype=np.uint)
        for (i, row) in enumerate(values):
            for (k, v) in group_by_value(row).iteritems():
                histogram[i][k] = v
        return histogram

#---------------------------------------------
# Random number generation
#---------------------------------------------
def random_index(n, size=None):
    '''Return an array of size 'size' containing distinct random integers between 0 and n-1. If
    size > n or size=None, n integers are returned (a random permutation of arange(n)).'''
    size = n if size is None else size
    return np.random.permutation(n)[:size]

def random_choice(a, n=None):
    '''Choose a random sample of the array (similar to the numpy 1.7.0 function choice, but we use
    numpy 1.6 right now.'''  
    return a[random_index(len(a), size=n)]

def rand_binary(n):
    '''Return a random binary sequence of size n.'''
    return np.array([np.random.randint(2) for _ in xrange(1, n + 1)])

def multinomial_elementwise(p, r=1):
    '''Return an m x r matrix x whose ith entry is a multinomial random variable in 0..n-1, distributed
    according to the probabilities in the ith row of the m x n row-stochastic matrix p. That is,
    p[i,j] is the probability that x[i,m] == j for every m=0..r-1. If r=1, x is an m-vector.'''
    (m, n) = p.shape
    # Create a 1-D vector of random[0,1] (t) and a placeholder for the integer multinomials (x)
    (x, t) = (np.zeros((r * m,), dtype=np.int), np.random.random((r * m,)))
    # Build cumulative distribution in each row
    p_cum = np.concatenate((np.zeros((r * m, 1)), np.tile(np.cumsum(p, axis=1), (r, 1))), axis=1)    
    for i in xrange(n):
        x[np.where((p_cum[:, i] <= t) & (t < p_cum[:, i] + p_cum[:, i + 1]))] = i
    return x if r == 1 else x.reshape((m, r), order='F')

def random_row_stochastic(size):
    '''Return a row-stochastic matrix of size size=(m,n).'''
    return scale_row_sums(np.random.random(size))

def scale_row_sums(a):
    '''Scale a matrix's row sums to 1.'''
    return a / a.sum(axis=1)[:, np.newaxis]

#---------------------------------------------
# Vector and Matrix Norms
#---------------------------------------------
def norm_frobenius_scaled(a):
    '''Return the scaled Frobenius norm of the array a.'''
    return np.linalg.linalg.norm(a, 'fro') / np.sqrt(a.size)

#---------------------------------------------
# Scaling
#---------------------------------------------
'''Functions to scale the a positive measure (e.g., POO phase measure) to [-1,1].'''
POSITIVE_TO_0_1 = lambda x: x / (x + 1.)
_t = lambda x: np.sign(x) * POSITIVE_TO_0_1(np.abs(x))
POSITIVE_TO_MINUS_1_TO_1 = lambda x: _t(np.log(x))

#---------------------------------------
# Mathematical & Statistical functions
#---------------------------------------
iabs = lambda x: x if x >= 0 else -x

def log_factorial(N):
    '''Return an array of log(n) values for n = 0..N.'''
    L = [0] * (N + 1)
    for n in xrange(2, N + 1): L[n] = L[n - 1] + np.log(n)
    return L

def log10_factorial(N):
    '''Return an array of log10(n) values for n = 0..N.'''
    L = [0] * (N + 1)
    for n in xrange(2, N + 1): L[n] = L[n - 1] + np.log10(n)
    return L

def binom():
    '''A generator of binomial coefficients. Iterate n returns [C(n,k) for k in xrange(n)], n = 0,1,2,... .'''
    b = [1L]
    yield list(b)
    for n in it.count(1):
        a = list(b)
        for i in xrange(n - 1): b[i] = a[i] + a[i + 1]
        b.insert(0, 1L)
        yield list(b)

def binom_mod(r):
    '''A generator of binomial coefficients modulo r >= 2.
    Iterate n returns [C(n,k) for k in xrange(n)], n = 0,1,2,... .'''
    b = [1L]
    yield list(b)
    for n in it.count(1):
        a = list(b)
        for i in xrange(n - 1): b[i] = (a[i] + a[i + 1]) % r
        b.insert(0, 1L)
        yield list(b)

def sum_mod(iterable, r):
    '''Summation modulo r of the iterates of iterable.'''
    return reduce(lambda x, y: (x + y) % r, iterable, 0)

def prod_mod(iterable, r):
    '''Product modulo r of the iterates of iterable.'''
    return reduce(lambda x, y: (x * y) % r, iterable)

def bernoulli(n, p, k):
    '''Probability of getting k successes in n Bernoulli trials with success probability p.'''
    if k < 0 or k > n: return 0
    if p == 0: return 1 if k == 0 else 0
    if p == 1: return 1 if k == n else 0
    return stats.binom(n, p).pmf(k)

def cumbin(n, p, m):
    '''Cumulative binomial disribution - probability of getting <= m successes in n Bernoulli trials
    with success probability p.'''
    return betainc(n - m, m + 1, 1 - p)
