'''
============================================================
A hashable decorator of a numpy array.

Using the wrapper class is simple enough:
>>> from numpy import arange
>>> a = arange(0, 1024)
>>> d = {}
>>> d[a] = 'foo'
TypeError: unhashable type: 'numpy.ndarray'
>>> b = hashable(a)
>>> d[b] = 'bar'
>>> d[b]
'bar'

Created on Septeber 14, 2012
@author: http://machineawakening.blogspot.com/2011/03/making-numpy-ndarrays-hashable.html
============================================================
'''
#from hashlib import sha1
from numpy import all, array #, uint8

class Hashable(object):
    r'''Hashable wrapper for ndarray objects.

        Instances of ndarray are not hashable, meaning they cannot be added to
        sets, nor used as keys in dictionaries. This is by design - ndarray
        objects are mutable, and therefore cannot reliably implement the
        __hash__() method.

        The hashable class allows a way around this limitation. It implements
        the required methods for hashable objects in terms of an encapsulated
        ndarray object. This can be either a copied instance (which is safer)
        or the original object (which requires the user to be careful enough
        not to modify it).
    '''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, wrapped, tight=False):
        r'''Creates a new hashable object encapsulating an ndarray.

            wrapped
                The wrapped ndarray.

            tight
                Optional. If True, a copy of the input ndaray is created.
                Defaults to False.
        '''
        self.__tight = tight
        if isinstance(wrapped, list):
            wrapped = array(wrapped)
        self.__wrapped = array(wrapped) if tight else wrapped
        #self.__hash = int(sha1(wrapped.view(uint8)).hexdigest(), 16) # Fast but unreliable
        #self.__hash = int(sha1(wrapped.astype(uint8)).hexdigest(), 16) # Reliable but slow
        self.__hash = hash(wrapped.tostring())
        #print 'Creating Hashable', self.__wrapped, self.__hash
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        r'''Return a textual representation of this object.'''
        #return 'Hashable[%s, %s]' % (repr(self.__wrapped.shape), repr(self.__wrapped))
        return 'H%s' % (repr(self.__wrapped),)

    def __eq__(self, other):
        return all(self.__wrapped == other.__wrapped)

    def __hash__(self):
        return self.__hash

    def unwrap(self):
        r'''Returns the encapsulated ndarray.

            If the wrapper is "tight", a copy of the encapsulated ndarray is
            returned. Otherwise, the encapsulated ndarray itself is returned.
        '''
        return array(self.__wrapped) if self.__tight else self.__wrapped