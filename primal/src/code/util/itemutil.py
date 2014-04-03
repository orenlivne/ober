'''
============================================================
Collection item extraction utilities.

Created on Jun 1, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np, itertools as it

#---------------------------------------------
# Methods
#---------------------------------------------
def ilen(it):
    '''Return the number of items in an iterator / generator expression.'''
    return sum(1 for _ in it)

def items(it, index):
    '''Return the items in an iterable by an index array.'''
    return (x for i, x in enumerate(it) if i in index)

def irange(it, start, stop, step=1):
    '''Return item range.'''
    return items(it, xrange(start, stop, step))

def linerange(it, start, stop, step=1):
    '''Return a line range in a file. Strip new lines.'''
    return (line.strip() for line in items(it, xrange(start, stop, step)))

def segmentrange(it, step):
    '''Return an Nx2 array of start and stop positions of N equidistant segments of size step into the
    item collection it.'''
    # Loop over items, save items that are alternatingly segment start and end points 
    (next_position, start, endpoints) = (0, True, [])
    for i, x in enumerate(it):
        if i == next_position:
            endpoints.append(x)
            next_position += ((step - 1) if start else 1)
            start = not start
    # x conveniently holds the last item; add it if the last segment incomplete
    if not start:
        endpoints.append(x)
    # Reshape into desired shape
    return np.array(endpoints).reshape((len(endpoints) / 2, 2))

'''Collect data into fixed-length chunks or blocks (possibly except the last).
grouper('ABCDEFG', 3) --> ABC DEF G'''
grouper = lambda iterable, n: (tuple(x for x in block if x) for block in it.izip_longest(*([iter(iterable)] * n))) 

def group_iterates(it, c, index=False):
    '''Groups ParamIterator iterator (it) results into blocks of size c each.
    If index=True, iterates are also enumerated within the block.'''
    buf = []
    for k, x in enumerate(it):
        buf.append((k, x) if index else x)
        if k % c == (c - 1):
            yield buf
            buf = []
    if buf:
        yield buf

def groupby_with_range(iterator):
    '''Group consecutive values and return their index range in the iterator. Ranges are
    [a,b], i.e., inclusive at both ends.
    
    Usage:
    
    itemutil.groupby_with_range('AAAABBBCCDAABBB') -->
    [('A', (0, 3)),
     ('B', (4, 6)),
     ('C', (7, 8)),
     ('D', (9, 9)),
     ('A', (10, 11)),
     ('B', (12, 14))]'''
    groups = zip(*((k, len(list(g))) for k, g in it.groupby(iterator)))
    if groups:
        (a, b) = tuple(groups)
        return zip(a, zip(np.cumsum(np.concatenate(([0], b[:-1]))), np.cumsum(b) - 1))
    else:
        return []

def ntuples(lst, n):
    '''Return periodic consecutive n-tuples.

    Usage:
    
    itemutil.ntuples(range(10), 3) -> 
    [(0, 1, 2),
     (1, 2, 3),
     (2, 3, 4),
     (3, 4, 5),
     (4, 5, 6),
     (5, 6, 7),
     (6, 7, 8),
     (7, 8, 9),
     (8, 9, 0),
     (9, 0, 1)]'''
    return it.izip(*[it.chain(it.islice(lst, i, None), it.islice(lst, None, i)) for i in xrange(n)])

def index_of_change(lst, output_value=None, output_first=False, output_last=False,
                    comparator=lambda x, y: x == y):
    '''Yield the indices of a list where its value changes. If output is specified, output
    the two iterates i, i+1 abutting the index of change i. If output_first=True, the first yielded
    entry is index=-1, item=the first item on the list, and prev=None.

    Examples:
    
    list(itemutil.index_of_change([0] + [1]*2 + [3]*3)) 
    --> [0, 2]
    
    list(itemutil.index_of_change([0] + [1]*2 + [3]*3, True)) 
    --> 
    [(0, 0, 1), (2, 1, 3)]
    
    list(itemutil.index_of_change([0] + [1]*2 + [3]*3, True, True, True))
    --> [(-1, None, 0), (0, 0, 1), (2, 1, 3), (4, 3, None)]   
    '''
    i = iter(lst)
    prev = item = i.next()
    if output_first:
        index = -1
        yield (index, None, item) if output_value else index
    for index, item in enumerate(i):
        if not comparator(prev, item):
            yield (index, prev, item) if output_value else index
        prev = item
    if output_last:
        yield (index, prev, None) if output_value else index
    
def binseq(k):
    '''All binary sequences of size k.'''
    for q in it.product('01', repeat=k):
        yield [int(x) for x in q]
    
def ranges(n, p):
    '''Break the integers 0..n-1 into p more-or-less equal ranges. This is an iterator over ranges.
    Usage:
    >>> ranges(10,3) = [(0, 3), (4, 6), (7, 9)]'''
    a = list(range_sizes(n, p))
    return zip(np.cumsum([0] + a[:-1]), np.cumsum(a) - 1)

def range_sizes(n, p):
    '''Break the integers 0..n-1 into p more-or-less equal ranges. This is an iterator over ranges.
    Usage:
    >>> range_sizes(10,3) = [4,3,3]'''
    q, r = n / p, n % p
    return it.chain(it.islice(it.repeat(q + 1), r), it.islice(it.repeat(q), p - r))

def range_iterator(iterable, block_sizes):
    '''Yield blocks of variable sizes specified by the list block_sizes.
    Usage:
    >>> list(range_iterator(xrange(10), [4,3,3]))
    [(0,1,2,3),(4,5,6),(7,8,9)]
    '''
    sz_iter = iter(block_sizes)
    b, count, sz = [], 0, sz_iter.next()
    for x in iterable:
        if count == sz:
            yield b
            b, count, sz = [], 0, sz_iter.next()
        print b, x
        b.append(x)
        count += 1
    if b:
        yield b
    
def block_iterator(iterable, block_size):
    '''Yield blocks of block_size iterates from the iterable iterable.
    Usage:
    >>> list(block_iterator(xrange(10), 3)) = [(1, 2), (3, 4, 5), (6, 7, 8), (9,)]'''
    return (tuple(y for y in x if y) for x in it.izip_longest(*([iter(iterable)] * block_size)))

#---------------------------------------------
# Private Methods
#---------------------------------------------
