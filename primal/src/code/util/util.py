'''
============================================================
Utility functions and classes.

Created on Jun 1, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import atexit, sys, tempfile, itertools as it, os, numpy as np, urllib2, StringIO, \
gzip, re, subprocess as sub, csv, cStringIO
from bst import BinarySearchTree
from scipy import ndimage
from scipy.spatial.distance import pdist

#---------------------------------------------
# Constants
#---------------------------------------------
'''Abnormal main program termination codes'''
EXIT_FAILURE = 1
EXIT_BAD_INPUT_ARGS = 2
EXIT_FILE_NOT_FOUND = 3

#---------------------------------------------
# Bit Operations
#---------------------------------------------
# For pretty-printing booleans
BOOL_STRING = {False: 'NO', True: 'YES'}

def hamming_distance(i, j): 
    '''Return the Hamming distance (number of different bits) between two uint8's i and j.'''
    i ^= j 
    i = ((i & 044) >> 2) + ((i & 0222) >> 1) + (i & 0111) 
    i = ((i & 070) >> 3) + (i & 0307)
    return int(i % 63)

def num_differences(x, y):
    '''Return the # differences between the between the BINARY sequences x and y.'''
    return pdist([x, y], 'hamming')[0] * len(x)

#---------------------------------------------
# Sets
#---------------------------------------------
def is_member(a, s):
    '''Check if all members of s are in the set a.'''
    return all(i in a for i in s)

def num_members(a, s):
    '''Output the number of members of s that are in the set a.'''
    return sum(i in a for i in s)

def has_at_least_n_members(a, s, n=None):
    '''Check if there are at least n members in the intersection(a,s). If n is None, checks whether
    all members of s are in a (i.e., n = len(s) effectively).'''
    if not n:
        return all(i in a for i in s)
    else:
        n = min(n, len(s))  # Can't expect to have more elements in intersection(a,s) than s's size
        count = 0
        for i in s:
            if i in a:
                count += 1
            if count == n:
                return True
        return False

def group_index(a):
    '''In a collection a, return the indices of the start of each group of identical
    consecutive elements.'''
    (i_prev, x_prev) = (-1, None)
    for (i, x) in enumerate(a):
        if i_prev < 0 or x_prev != x:
            (i_prev, x_prev) = (i, x)
            yield i

def remove_duplicates(haps):
    '''Given a set of key-value tuples, remove all tuples that share the same key.'''
    key_count = {}
    for (k, _) in haps:
        key_count.setdefault(k, 0)
        key_count[k] += 1
    return set((k, v) for (k, v) in haps if key_count[k] == 1) 

#    values.sort()
#    edges = (values[1:] != values[:-1]).nonzero()[0] - 1
#    idx = np.concatenate(([0], edges, [len(values)]))
#    index = np.empty(len(idx) - 1, dtype= 'u4,u2')
#    index['f0'] = values[idx[:-1]]
#    index['f1'] = np.diff(idx)
#    return index

def sort_with_index(a):
    '''Sort an array and return the index, i.e., return the sorted b and i such that b = a[i].'''
    i = sorted(range(len(a)), key=a.__getitem__)
    return (a[i], i)

#---------------------------------------------
# Dictionaries
#---------------------------------------------
class mdict(dict):
    '''Multivalued dictionary from 
    http://code.activestate.com/recipes/440502-a-dictionary-with-multiple-values-for-each-key/'''
    
    '''Add the given value to the list of values for this key'''
    def __setitem__(self, key, value): self.setdefault(key, []).append(value)

    @staticmethod
    def from_items(items):
        '''Initialize from a key-value tuple list.'''
        a = mdict()
        for k, v in items: a[k] = v
        return a
    
'''Merge two dictionaries. From http://stackoverflow.com/questions/38987/how-can-i-merge-union-two-python-dictionaries-in-a-single-expression'''
dictmerge = lambda x, y: dict(x.items() + y.items())

'''Convert a key-value-pair list a=[ (1, 'A'), (1, 'B'), (2, 'C') ] into a multi-valued dictionary,
so that all values with the same key would be aggregated into a set:
{ 1: set(['A', 'B']), 2:(set['C'],) }. The dictionary is sorted by keys.'''
to_set_dict = lambda a: dict((key, set(v for (_, v) in pairs)) for key, pairs in it.groupby(sorted(a), lambda pair: pair[0]))

def occur_dict(items):
    '''A dictionary of item occurrences in the iterable items.'''
    d = {}
    for i in items:
        if i in d: d[i] = d[i] + 1
        else: d[i] = 1
    return d

'''Invert a dictionary.'''
dict_invert = lambda d: dict(zip(d.itervalues(), d.iterkeys()))

def dict_to_array(d):
    '''Convert a dictionary to numpy array.'''
    a = np.empty(len(d), dtype={'names': ('k', 'v'), 'formats': ('u2', 'u2')})
    a['k'] = d.keys()
    a['v'] = d.values()
    return a

'''Return the union of a list of sets.'''
union_all = lambda *sets: set().union(*sets)

#---------------------------------------------
# Binary Search Trees and Nearest Neighbors
#---------------------------------------------
def alternating_order(n):
    '''Return a numpy array with a permutation of 0..n-1 that starts in the middle and 
    alternates between left and right steps. Useful for efficient BSD insertion order of
    an ordered list. But much better is optimal_insertion_order().'''
    result = np.empty(n, dtype=int)
    if np.mod(n, 2) == 0:
        # Even n
        result[0:n:2] = np.arange(n / 2, n)
        result[1:n:2] = np.arange(n / 2 - 1, -1, -1) 
    else:
        # Odd n
        result[0:n:2] = np.arange((n - 1) / 2, n)
        result[1:n:2] = np.arange((n - 3) / 2, -1, -1) 
    return result

def optimal_insertion_order(n):
    '''Best insertion order of a sorted list into a BST. Forces all left subtrees to be full.'''
    depth = int(np.floor(np.log2(n))) + 1
    order = np.empty(n, dtype=int)
    start = 2 ** (depth - 1)
    step = 2 * start
    counter = 0
    for _ in xrange(0, depth):
        a = np.arange(start - 1, n, step)
        added = np.size(a)
        counter_new = counter + added
        order[counter:counter_new] = a
        counter = counter_new
        step = step / 2
        start = start / 2 
    return order

def sequence_to_tree(values, sort_key=None):
    '''Convert a sorted numpy array 'values' into a tree using the optimal insertion order.'''
    return BinarySearchTree(sort_key=sort_key, values=values[optimal_insertion_order(np.size(values))])

def list_index_tree(locations):
    '''Convert a sorted numpy array 'values' into a BST tree.'''
    last = len(locations) - 1
    index = np.arange(0, last + 1)
    return sequence_to_tree(index, sort_key=dict(zip(index, locations)).__getitem__)

def nearest_neighbor(x, locations):
    '''Return the nearest location to x in the list locations.'''
    return nearest_neighbor_in_list_tree(x, locations, list_index_tree(locations))

def nearest_neighbor_multiple(x, locations):
    '''Yield the nearest location to each item of x in the list locations.'''
    return nearest_neighbor_in_list_tree_multiple(x, locations, list_index_tree(locations))

def nearest_neighbor_in_list_tree(x, locations, b):
    '''Yield the nearest location to x in the list locations. b is a prepared
    BST of the list.'''
    return __nearest_neighbor_in_list_tree_helper(x, locations, b, len(locations) - 1)

def nearest_neighbor_in_list_tree_multiple(x, locations, b):
    '''Yield the nearest location to each item of x in the list locations. b is a prepared
    BST of the list.'''
    last = len(locations) - 1
    for y in x:
        yield __nearest_neighbor_in_list_tree_helper(y, locations, b, last)

def __nearest_neighbor_in_list_tree_helper(x, locations, b, last):
    left_nbhr = b.find_largest_le(x)
    if left_nbhr is None:
        return 0
    elif left_nbhr == last:
        return last
    else:
        right_nbhr = left_nbhr + 1
        (left, right) = (locations[left_nbhr], locations[right_nbhr])
        return left_nbhr if x - left < right - x else right_nbhr

def nearest_neighbor_in_sorted_arrays(a, b):
    '''For each element of the sorted array b, emit the index of the closest element of a.''' 
    i, n = 0, len(a)
    for x in b:
        while i < n and a[i] < x:
            i += 1
        yield i - 1 if (i > 0) and ((i == n) or (x - a[i - 1] <= a[i] - x)) else i

#---------------------------------------------
# Numpy Grids and N-dimensional Arrays
#---------------------------------------------
def brange(start, stop, step=1, endpoint=False, dtype=None):
    '''Like numpy's arange, but adds the option of retaining the last endpoint if endpoint is True
    and it is not already included in the arange.'''
    a = np.arange(start, stop, step=step, dtype=dtype)
    if endpoint and a[-1] != stop:
        a = np.append(a, stop) 
    return a
 
def flattened_meshgrid(x, y):
    '''Return flattened meshgrid arrays, like a Kronecker product of x and y.'''
    return [a.flatten() for a in np.meshgrid(x, y)]

def lexsort_by_columns(a):
    '''Indirectly-lexicographically-sort a numpy matrix a by column 0, then by column 1, ..., 
    column a.shape[1]-1.'''
    return np.lexsort([x for x in np.flipud(a.transpose())])

def nans(shape, dtype=float):
    '''Return an array filled with NaN values.'''
    a = np.empty(shape, dtype=dtype)
    a.fill(np.nan)
    return a

def sum_by_group(values, groups):
    '''Sum a numpy array ''values'' by group. Group indices are specified by the ''groups'' array.
    @see http://stackoverflow.com/questions/4373631/sum-array-by-number-in-numpy''' 
    order = np.argsort(groups)
    groups = groups[order]
    values = values[order]
    values = values.cumsum()
    index = np.ones(len(groups), 'bool')
    index[:-1] = groups[1:] != groups[:-1]
    values = values[index]
    groups = groups[index]
    values[1:] = values[1:] - values[:-1]
    return values, groups

#---------------------------------------------
# Numpy - Misc
#---------------------------------------------
def set_printoptions(options):
    '''Set print options from a struct obtained by get_printoptions().'''
    np.set_printoptions(edgeitems=options['edgeitems'],
                        infstr=options['infstr'],
                        linewidth=options['linewidth'],
                        nanstr=options['nanstr'],
                        precision=options['precision'],
                        suppress=options['suppress'],
                        threshold=options['threshold'])

#---------------------------------------------
# Digital Filters
#---------------------------------------------
def max_except_center_filter(a, size):
    '''Applies the following digital filter to a {-1,0,1}-valued array a: If a is not 1, return
    0; otherwise return the maximum over the window [-(size-1)/2, (size-1)/2] with the central
    element excluded.'''
    middle = (size - 1) / 2
    footprint = range(size)
    footprint.remove(middle)        
    return ndimage.generic_filter(a, __max_except_center, size, extra_arguments=(footprint, middle),
                                  mode='constant', cval=0)

def max_except_center_filter1d(a, size, axis= -1):
    '''Applies the following digital filter to a {-1,0,1}-valued array a: If a is not 1, return
    0; otherwise return the maximum over the window [-(size-1)/2, (size-1)/2] with the central
    element excluded.'''
    return ndimage.generic_filter1d(a, __max_except_center_1d, size, extra_arguments=(size,),
                                    axis=axis, mode='constant', cval=0)

def __max_except_center(window, footprint, middle):
    '''Assuming symmetric filter footprint.'''
    return 0 if window[middle] != 1 else max(window[footprint])

def __max_except_center_1d(input_line, output_line, size):
    '''Assuming symmetric filter footprint.'''
    output_line[:] = max_except_center_filter(input_line[(size + 1) / 2:(len(input_line) - (size - 3) / 2)], size)

#---------------------------------------------
# I/O
#---------------------------------------------
'''Return the extension of a file string.'''
file_ext = lambda filename: os.path.splitext(filename)[1][1:]

def mkdir_if_not_exists(path):
    '''Create a directory if it does not exist yet.'''
    # print 'Making directory', path
    if path and not os.path.exists(path): os.makedirs(path)

def temp_filename(prefix=None, suffix='tmp', dir=None, text=False, removeOnExit=True):  # @ReservedAssignment
    ''''Returns a temporary filename that, like mkstemp(3), will be secure in
    its creation.  The file will be closed immediately after it's created, so
    you are expected to open it afterwards to do what you wish.  The file
    will be removed on exit unless you pass removeOnExit=False.  (You'd think
    that amongst the myriad of methods in the tempfile module, there'd be
    something like this, right?  Nope.)'''
    if prefix is None: prefix = "%s_%d_" % (os.path.basename(sys.argv[0]), os.getpid())
    fileHandle, path = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=dir, text=text)
    os.close(fileHandle)

    def removeFile(path):
        os.remove(path)
        # logging.debug('temporaryFilename: rm -f %s' % path)

    if removeOnExit: atexit.register(removeFile, path)
    return path

def check_file_is_readable(path):
    '''Check if file is readable; if not, throw an exception.'''
    try:
        with open(path, 'rb'): pass
    except IOError:
        raise IOError('Could not open file %s for reading' % (path,))
       
def check_file_is_writable(path):
    '''Check if file is readable; if not, throw an exception. Note: the file will be removed
    if this method is called and no exception is raised.'''
    try:
        with open(path, 'wb'): pass
        os.remove(path)
    except IOError:
        raise IOError('Could not open file %s for writing' % (path,))

def block_readlines(path, buffer_size):
    '''An iterator that reads a block of lines from a file at a time. buffer_size is the buffer
    size, not the number of lines in a block.'''
    with open(path, 'rb') as f:
        lines = f.readlines(buffer_size)
        while lines:
            yield [line.rstrip('\r\n').rstrip('\n') for line in lines]
            lines = f.readlines(buffer_size)
    
def module_path():
    '''Return the path to the module of the currently executing python file. Unlike the executing
    file path, which is indeterminate, the module directory is fixed.
    
    This of course means that if you call this method from another file, it will return the module
    of THIS file, not the calling file.
    
    Code taken from http://stackoverflow.com/questions/2632199/how-do-i-get-the-path-of-the-current-executed-file-in-python'''
    def __we_are_frozen():
        # All of the modules are built-in to the interpreter, e.g., by py2exe
        return hasattr(sys, "frozen")
    encoding = sys.getfilesystemencoding()
    return os.path.dirname(unicode(sys.executable, encoding)) if __we_are_frozen() \
        else os.path.dirname(unicode(__file__, encoding))

def run_command(cmd, die_on_exception=True, verbose=False):
    '''Run command in a sub-shell. Return the command's return code.'''
    if verbose:
        print 'Running', cmd
    # subprocess.call(cmd)
    p = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
    output, errors = p.communicate()
    if p.returncode != 0:
        if verbose:
            print 'Output', output
            print 'Errors', errors
        if die_on_exception:
            raise ValueError('Command failed:' , cmd)
    return p.returncode, output, errors

'''Return a generator expression of the lines in the delimiter-delimited text file f (a file
handle or file name), skipping the first num_header_lines.'''
read_file = lambda f, num_header_lines = 0, delimiter = '\t': \
    it.islice((line for line in csv.reader(f if isinstance(f, file) else open(f, 'rb'), delimiter=delimiter, skipinitialspace=True) if line),
              num_header_lines, None)

def write_buffered_lines(lines, out, buf_size=10000, delimiter='\t'):
    '''Write a list of lines of items to the output stream out, delimited by delimiter.
    Buffer the output.'''
    # Read # samples in the first relevant line
    i, buf = 0, cStringIO.StringIO()
    for line in lines:
        # Write line to buffer
        i += 1
        buf.write(delimiter.join(line) + '\n')

        # Flush buffer if overflows
        if i == buf_size:
            # Write buffer to output stream; clear buffer
            out.write(buf.getvalue())
            buf.close()
            buf = cStringIO.StringIO()
            i = 0
        
    # Write the last buffer if it is non-empty to the output stream
    if i: out.write(buf.getvalue())
    buf.close()

#---------------------------------------------
# I/O - Networking
#---------------------------------------------
__URL_OPENER = urllib2.build_opener()

def url_stream(url):
    '''Download and uncompress a gzipped web URL; it is returned as a file stream.'''  
    request = urllib2.Request(url)
    request.add_header('Accept-encoding', 'gzip')
    return __URL_OPENER.open(request)

def url_gzip_stream(url):
    '''Download and uncompress a gzipped web URL; it is returned as a file stream.'''  
    request = urllib2.Request(url)
    request.add_header('Accept-encoding', 'gzip')
    f = __URL_OPENER.open(request)
    compresseddata = f.read()
    compressedstream = StringIO.StringIO(compresseddata)
    return gzip.GzipFile(fileobj=compressedstream)

__URL_PATTERN = re.compile('(\w+)://(.+)')

def open_resource(url):
    '''Open URL/local file as an input stream.'''
    if url.endswith('.gz'):
        print 'Reading zipped URL %s ...' % (url,) 
        return url_gzip_stream(url)
    else:
        m = __URL_PATTERN.match(url)
        protocol = m.group(1)
        local_file = not m or (protocol == 'file')
        if local_file:
            file_url = m.group(2) if m and local_file else url 
            print 'Reading local file %s ...' % (file_url,) 
            return open(file_url, 'rb')
#        elif protocol == 'ssh':
#            print 'Reading remote SSH URL %s ...' % (url,) 
#            return open(url_stream(url), 'rb')            
        else:
            print 'Reading remote URL %s ...' % (url,) 
            return open(url_stream(url), 'rb')

#---------------------------------------------
# Printouts
#---------------------------------------------
def cprint(a, fmt, prefix='[', postfix=']'):
    '''Print a collection in a printf-style format ''fmt''.'''
    return prefix + ' '.join(fmt % x for x in a) + postfix

#---------------------------------------------
# Parallel Processing
#---------------------------------------------
def write_with_lock(s, lock):
    '''Write the string ''s'' to stdout; synchronize if a parallel pool lock is specified.'''
    if lock:
        lock.acquire()
        sys.stdout.write(s)
        sys.stdout.flush()
        lock.release()
    else: sys.stdout.write(s)

####################################################################################
class Struct(object):
    '''A generic python struct class. Taken from http://norvig.com/python-iaq.html'''

    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, **entries): self.__dict__.update(entries)

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self): return repr(self.__dict__) 

    '''Override existing struct values with new values.'''
    def update(self, **entries): self.__dict__.update(entries)

    '''Override existing struct values with new values from another Struct instance.'''
    def update_from_struct(self, other): self.__dict__.update(other.__dict__)
