'''
============================================================
Fibonacci number calculations. 
============================================================
'''
import timeit
from exp import bits

def fib1(n):
    '''Recursive. Exponential runtime.'''
    # print n
    if n == 0: return 0
    elif n == 1: return 1
    else: return fib1(n - 1) + fib1(n - 2)

def fib2(n, f0=0, f1=1):
    '''Iterative. Uses an array for repeated subproblems'''
    for _ in xrange(2, n + 1):
        f0, f1 = f1, f0 + f1
    return f1

def matrix_mult(x, y):
    return (x[0] * y[0] + x[1] * y[2], x[0] * y[1] + x[1] * y[3],
            x[2] * y[0] + x[3] * y[2], x[2] * y[1] + x[3] * y[3])

def matrix_pow(x, n):
    '''D & C, iterative.'''
    y = (1, 0, 0, 1)
    for r in reversed(list(bits(n))):
        y = matrix_mult(y, y)
        if r: y = matrix_mult(y, x)
    return y

def fib3(n):
    '''Divide and conquer using 2x2 matrix exponentiation.'''
    return matrix_pow((1, 1, 1, 0), n)[1]

def time_fib(impl, k_max=20):
    n = 2
    for k in xrange(1, k_max + 1):
        print n, timeit.timeit('%s(%d)' % (impl, n,), 'from __main__ import %s' % (impl,), number=5)
        n *= 2
    
if __name__ == "__main__":
    print '--- fib2 ---'
    time_fib('fib2')
    print '--- fib3 ---'
    time_fib('fib3')
    #for n in xrange(10):
        #print n, fib2(n), fib3(n)
    # time_fib('fib1')
