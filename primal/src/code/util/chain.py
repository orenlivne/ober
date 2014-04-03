'''
============================================================
Phasing chain framework.

The phasing algorithm is a chain comprising of filters.
In our context, a phasing stage (=filter) runs on a
haplotype (=the request) and passes it on to the next 
filters.

Created on July 5, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
#import numpy as np
#from pedigree.Person import Person

####################################################################################
class Filter(object):
    '''A filter abstraction.'''
    def __init__(self, next_filter=None, name=None, handle=None, debug=False):
        '''Initialize a filter. Specify its next filter.'''
        self.next_filter = next_filter
        self.name = name if name is not None else self.__class__.__name__
        self._debug = debug
        if handle is not None: self.handle = lambda y: handle(self, y)

    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def __repr__(self):
        return 'Filter[%s]' % (self.name, ) 
     
    def handle(self, request):
        '''Run the filter on a request. Return a 'handled' code: iff True,
        the request will not be passed on to self.next_filter.'''
        raise ValueError('Sub-classes must implement this method') 
    
    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, debug):
        self._debug = debug
                 
####################################################################################
class FilterChain(Filter):
    '''A generic filter chain.'''
    # TODO: make debug=True act like AOP to automatically decorate to filter decorator delegates
    # (e.g., SnpTrioPhaser.processor) 
    def __init__(self, chain, next_filter=None, debug=False, name=None):
        '''Initialize a filter chain from a list of filters. If the debug flag is on, prints
        which filters are invoked during request handling.'''
        self.chain = chain
        # Decorate if necessary
        if debug:
            self.chain = [PrintDecorator(f) for f in self.chain]

        Filter.__init__(self, next_filter, name=name, debug=debug)
         
        # Connect each consecutive filters  
        prevFilter = None 
        for f in self.chain:
            if prevFilter:
                f.next = prevFilter
            prevFilter = f
        
        self.debug = debug
        
    def handle(self, request):
        '''Run the filter chain on a request.'''
        for f in self.chain:
            if f.handle(request):
                # Premature chain termination when a chain link handles the request
                return True
        # Entire chain is done but may be embedded within a bigger chain
        return False

    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, debug):
        self._debug = debug
        # Propagate debug mode to child filters
        for f in self.chain:
            f.debug = debug

####################################################################################
class PrintFilter(Filter):
    '''An example filter that prints its name and the request.'''
    def __init__(self, name, next_filter=None):
        '''Initialize a name-printing filter.'''
        Filter.__init__(self, next_filter)
        self.name = name
     
    def handle(self, request):
        '''Print some info.'''
        print '[' + self.name + ']: ' + repr(request) 

####################################################################################
class PrintDecorator(Filter):
    '''A filter decorator that prints a name upon request handling.'''
    def __init__(self, delegate, next_filter=None, fmt='[%-22s]'):
        '''Initialize a decorator of the delegate Filter. fmt specifies the printout format string that
        accepts one string argument for the filter's name.'''
        self.delegate = delegate
        self.fmt = fmt
        Filter.__init__(self, next_filter, name=delegate.name)
     
    def handle(self, request):
        '''Deocrate request handling without printouts.'''
        print (self.fmt + ' START') % (self.name,)
        handled = self.delegate.handle(request)
        print (self.fmt + ' END') % (self.name,)
        return handled 

    @property
    def debug(self):
        return self._debug

    @debug.setter
    def debug(self, debug):
        self._debug = debug
        self.delegate.debug = debug
