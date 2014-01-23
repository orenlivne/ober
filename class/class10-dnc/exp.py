'''
============================================================
Exponentiation algorithms.
============================================================
'''
from numpy.ma.testutils import assert_equal

def pow1(x, n):
    return 1 if n == 0 else x * pow1(x, n - 1)

def pow2(x, n):
    y = 1
    for i in xrange(n): y *= x
    return y

def pow2a(x, n):
    return reduce(lambda y, _: y * x, xrange(n), 1)

def pow3(x, n):
    if n == 0: return 1
    else:
        n, r = divmod(n, 2)
        y = pow3(x, n)
        y = y * y
        if r: y *= x
        return y

def bits(n):
    '''Return the bits in the binary representation of n, left-to-right.'''
    while n:
        n, r = divmod(n, 2)
        yield r

def pow4(x, n):
    '''D & C, iterative.'''
    y = 1
    for r in reversed(list(bits(n))):
        y *= y
        if r: y *= x
    return y

def test_pow():
    x = 123
    n = 16
    expected = x ** n
    assert_equal(pow1(x, n), expected, 'Wrong decrease & conquer recursive implementation')
    assert_equal(pow2(x, n), expected, 'Wrong decrease & conquer iterative implementation')
    assert_equal(pow2a(x, n), expected, 'Wrong decrease & conquer reduce implementation')
    assert_equal(pow3(x, n), expected, 'Wrong D & C, recursive implementation')
    assert_equal(pow4(x, n), expected, 'Wrong D & C, iterative implementation')
    
if __name__ == "__main__":
    test_pow()
