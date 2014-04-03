'''
============================================================
Created on Jun 1, 2012

@author: Oren Livne <livne@uchicago.edu>
============================================================
'''

class Angle(object):
    def __init__(self, rad):
        self._rad = rad

    @property
    def rad(self):
        '''The angle in radians'''
        return self._rad
        
    @rad.setter
    def rad(self, value):
        self._rad = float(value)

    @rad.deleter
    def rad(self):
        del self._rad
  
####################################################################################
class Clazz2(object):
    '''Example, copied verbatim from http://docs.python.org/library/functions.html#property'''
    def __init__(self):
        self._x = 10

    def getx(self):
        return self._x
    def setx(self, value):
        self._x = value
    def delx(self):
        del self._x
    x = property(getx, setx, delx, 'I''m the ''x'' property.')
    
class Parrot(object):
    '''Example, copied verbatim from http://docs.python.org/library/functions.html#property'''
    def __init__(self):
        self._voltage = 100000

    @property
    def voltage(self):
        '''Get the current voltage.'''
        return self._voltage
    
####################################################################################
class Clazz1(object):
    '''Illustrates a property with a getter and setter.'''
    
    def __init__(self):
        #print 'Creating instance of Clazz1'
        self._x = None
        self.x = 200
        self._voltage = 100000

    @property
    def voltage(self):
        '''Get the current voltage.'''
        return self._voltage

    @property
    def x(self):
        '''I'm the 'x' property.'''
        return self._x

    @x.setter
    def x(self, value):
        self._x = value

    @x.deleter
    def x(self):
        del self._x

####################################################################################
class Seq:
    '''An iterator class that spits consecutive perfect squares.'''
    def __init__(self):
        self.x = 0

    def next(self): #@ReservedAssignment
        self.x += 1
        return self.x*self.x

    def __iter__(self):
        return self

####################################################################################
class Slicer(object):
    '''Slices x s.t. the index of the first dimension is a fixed scalar.'''
    def __init__(self, x, index):
        self.x = x
        self.index = index
    
    def __getitem__(self, *args):
        return self.x.__getitem__((self.index,) + args[0])
         
        #return np.logical_and(g[:,:,PATERNAL] != MISSING, g[:,:,PATERNAL] == g[:,:,MATERNAL])
    
####################################################################################
if __name__ == '__main__':
    '''Usage examples.'''
    angle = Angle(10)
    print angle.rad
    
    c = Clazz2()
    print 'x', c.x
    
    p = Parrot()
    print 'voltage', p.voltage
    
    c = Clazz1()
    print 'voltage', c.voltage, 'x', c.x
    c.x = 10
    print 'voltage', c.voltage, 'x', c.x
    assert c.voltage == 100000, 'Property not set properly'

    s = Seq()
    n = 0
    for i in s:
        print i
        n += 1
        if n > 10:
            break
