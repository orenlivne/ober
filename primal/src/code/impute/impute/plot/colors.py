'''
============================================================
Miscellaneous plot functions.

Created on Jan 14, 2013
@author: http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors 
original version for Python 3.3
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import colorsys, itertools, numpy as np, matplotlib as P
from fractions import Fraction

def zenos_dichotomy():
    '''
    http://en.wikipedia.org/wiki/1/2_%2B_1/4_%2B_1/8_%2B_1/16_%2B_%C2%B7_%C2%B7_%C2%B7
    '''
    for k in itertools.count():
        yield Fraction(1, 2 ** k)

def getfracs():
    '''
    [Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4), Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8), Fraction(1, 16), Fraction(3, 16), ...]
    [0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625, 0.1875, ...]
    '''
    yield 0
    for k in zenos_dichotomy():
        i = k.denominator  # [1,2,4,8,16,...]
        for j in xrange(1, i, 2): yield Fraction(j, i)

'''Can be used for the v in hsv to map linear values 0..1 to something that looks equidistant.'''
bias = lambda x: (np.sqrt(x / 3) / Fraction(2, 3) + Fraction(1, 3)) / Fraction(6, 5)

def genhsv(h):
    for s in [Fraction(6, 10)]:  # optionally use range
        for v in [Fraction(8, 10), Fraction(5, 10)]:  # could use range too
            yield (h, s, v)  # use bias for v here if you use range

genrgb = lambda x: colorsys.hsv_to_rgb(*x)
flatten = itertools.chain.from_iterable
gethsvs = lambda : flatten(itertools.imap(genhsv, getfracs()))
getrgbs = lambda : itertools.imap(genrgb, gethsvs())

fixed = ['r', 'g', 'b', 'c', 'm', 'y']
num_fixed = len(fixed)
fixed_rgb = lambda x: map(P.colors.colorConverter.to_rgb, x)
get_colors = lambda: itertools.chain(fixed_rgb(fixed), getrgbs())

def genhtml(x):
    uint8tuple = itertools.imap(lambda y: int(y * 255), x)
    return 'rgb({},{},{})'.format(*uint8tuple)

gethtmlcolors = lambda: itertools.imap(genhtml, get_colors())

if __name__ == '__main__':
    print '\n'.join(map(lambda x: '<div style="color:%s;">%s</div>"' % (x, x), itertools.islice(gethtmlcolors(), 100)))
