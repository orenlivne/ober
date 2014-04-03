'''
============================================================
Test singleton class decoration.

Created on April 2, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest

def singleton(cls):
    instances = {}
    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance

# def singleton(cls):
#     obj = cls()
#     cls.__new__ = staticmethod(lambda cls: obj)  # Always return the same object
#     return cls

# Can't have sub-classes

@singleton
class Duck(object):
    pass

# @singleton
# class SpecialDuck(Duck):
#     pass
    
class TestSingleton(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------   
    
    def test_duck_is_singleton(self):
        '''Test whether an method object argument is passed by reference.'''
        assert Duck() is Duck(), 'Singleton doesn''t work'
#         
#     def test_subclass_is_singleton(self):
#         '''Test whether an method object argument is passed by reference.'''
#         print id(Duck())
#         print id(SpecialDuck())
#         print id(SpecialDuck())
#         assert SpecialDuck() is SpecialDuck(), 'Singleton doesn''t work in subclass'
