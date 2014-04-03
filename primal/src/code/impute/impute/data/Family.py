'''
============================================================
A nuclear family of samples.

Created on July 22, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from impute.data import constants

class Family(object):
    #---------------------------------------------
    # Constructors
    #---------------------------------------------
    def __init__(self, parents, children=None):
        '''Construct a family from parents and a list of children.'''
        # Input validation.
        self._parents = parents
        self._children = sorted(set(children)) if children is not None else set([])
        
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    
    def __key(self):
        '''Uniquely-identifying key of this family: parent tuple.'''
        return (self.parents[constants.PATERNAL], self.parents[constants.MATERNAL])

    def __eq__(self, y):
        '''Equality of families.'''
        return self.__key() == y.__key()

    def __ne__(self, y):
        '''Inequality of families.'''
        return self.__key() != y.__key()

    def __hash__(self):
        '''Family hash code.'''
        return hash(self.__key())
        
    def __repr__(self):
        return 'Family[parents=(%d,%d), children=(%s)]' \
            % (self.parents[constants.PATERNAL], self.parents[constants.MATERNAL], ','.join(repr(x) for x in self.children), ) 
    
    def parents_dict(self):
        '''Return a dictionary of parent type (PATERNAL/MATERNAL) to parent node ID of node i.'''
        return dict(zip(constants.ALLELES, self._parents))
        
    def add(self, child):
        '''Add a child.'''  
        self._children.add(child) 

    def remove(self, child):
        '''Remove a child.'''  
        self._children.remove(child) 

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def parents(self):
        '''Parents' index tuple.'''  
        return self._parents

    @property
    def father(self):
        '''Father's index tuple.'''  
        return self._parents[constants.PATERNAL]
    
    @property
    def mother(self):
        '''Father's index tuple.'''  
        return self._parents[constants.MATERNAL]
    
    @property
    def children(self):
        '''Children's index set.'''  
        return self._children

    @property
    def num_children(self):
        '''# Children.'''  
        return len(self._children)

    @property
    def member_set(self):
        '''All family members: parents + children (set).'''  
        return  set(self._parents) | self._children

    @property
    def member_list(self):
        '''All family members: father, mother, children (list in this order).'''  
        return  list(self._parents) + list(self._children)
    
    @property
    def children_list(self):
        '''Children's list.'''  
        return list(self._children)
    
    @property
    def children_array(self):
        '''Children's list.'''  
        return np.array(self.children_list)
