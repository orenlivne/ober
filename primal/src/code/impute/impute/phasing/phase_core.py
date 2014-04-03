'''
============================================================
Phasing algorithm - resuable filter types.
 
Created on July 26, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import time, util
from chain import Filter, FilterChain
from impute.tools.param import PhaseParam
from impute.tools import pedigree_tools as pt

####################################################################################
class __PhaserChain(FilterChain):
    '''A base class for phasing processing chains.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, chain, name, debug, print_times, next_filter=None):
        '''Initialize a filter chain from a list of filters. Decorate them if debugging is on.
        Note: do not call outside this class. Use new_phaser_chain() instead.'''
        FilterChain.__init__(self, [PhaseDecorator(x) for x in chain] if print_times else chain,
                             name=name, debug=debug, next_filter=next_filter)
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def run(self, problem, params=None):
        '''Run the phasing processing chain. Adapts the generic Filter interface to include
        both a Problem and PhasingParam inputs. If params=None, using default PhaseParam values.''' 
        '''A template method that delegates to runner(), which accepts two input parameters.''' 
        return self.handle(util.Struct(problem=problem, params=params if params else PhaseParam()))

def new_phaser_chain(chain, name='Phasing', debug=False, print_times=False, next_filter=None):
    '''A factory method of phasing chains. Call this instead of __init__.'''
    chain = __PhaserChain(chain, name, debug=debug, print_times=print_times, next_filter=next_filter)
    return PhaseDecorator(chain) if print_times else chain

####################################################################################
class FamilyPhaser(Filter):
    '''A phaser template that sequentially processes families.'''
    
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, family_handle, next_filter=None, name=None,
                 genotyped_families=True, min_children=0):
        '''Initialize a filter that uses the filter processor to process SNP+trio combinations.
        Use single_family in debugging to process only the first encountered family.'''
        Filter.__init__(self, name=name, next_filter=next_filter)
        self.family_handle = family_handle
        self.genotyped_families = genotyped_families
        self.min_children = min_children
    
    #---------------------------------------------
    # Impl: Filter
    #---------------------------------------------
    def handle(self, request):
        '''Loop over all families and handle each one.'''
        problem, params = request.problem, request.params
        families = problem.find_families_by_member(params.single_member, genotyped=False, min_children=self.min_children) \
        if params.single_member else pt.selected_families(problem, params, genotyped=self.genotyped_families, min_children=self.min_children)
        for i, family in enumerate(families):
            # if self.debug and np.mod(snp, problem.num_snps/10):
            if params.debug:
                print 'Processing %d/%d %s' % ((i + 1), len(families), repr(family))
            # Note: ignoring processor handled (return) code
            self.family_handle(self, request, family)
            
            # Exit after one family -- debugging
            if params.single_family:
                return True
        return False

####################################################################################
class PhaseDecorator(Filter):
    '''A filter decorator of the main phasing steps.'''
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, delegate, next_filter=None):
        '''Initialize a decorator of the delegate Filter. fmt specifies the printout format string that
        accepts one string argument for the filter's name.'''
        Filter.__init__(self, next_filter=next_filter, name=delegate.name)
        self.delegate = delegate
     
    #---------------------------------------------
    # Impl: Filter
    #---------------------------------------------
    def run(self, problem, params=None):
        '''Run the phasing processing chain. Adapts the generic Filter interface to include
        both a Problem and PhasingParam inputs. If params=None, using default PhaseParam values.''' 
        '''A template method that delegates to runner(), which accepts two input parameters.''' 
        return self.handle(util.Struct(problem=problem, params=params if params else PhaseParam()))

    def handle(self, request):
        '''Deocrate problem handling without printouts.'''
        problem = request.problem
        h = problem.haplotype
        filled_before = 100.*h.fill_fraction()
        errors_before = 100.*problem.num_errors / h.num_data
        start = time.time()
        handled = self.delegate.handle(request)
        time_elapsed = time.time() - start
        filled_after = 100.*h.fill_fraction()
        errors_after = 100.*problem.num_errors / h.num_data
        print('[%-22s] Filled %-8d %5.2f%% (+%5.2f%%)  Errors %-5d %5.3f%% (+%5.3f%%)  Time %-5.1fs (%.1e s/data)' % \
              (self.name, h.num_filled, filled_after, filled_after - filled_before,
               problem.num_errors, errors_after, errors_after - errors_before,
               time_elapsed, 1.*time_elapsed / h.num_data,))
        return handled 

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def debug(self):
        return self.delegate.debug

    @debug.setter
    def debug(self, debug):
        self.delegate.debug = debug
