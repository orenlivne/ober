'''
============================================================
Test that gs is on the path. Prerequisite for using PIL 
image conversion functions.

Created on August 27, 2012
@author: Oren Livne <livne@uchicago.edu>
@see http://code.google.com/p/elaphe/issues/detail?id=7
============================================================
'''
from StringIO import StringIO
from PIL.EpsImagePlugin import EpsImageFile # May be an import error in Eclipse due to PIL's statically-linked imports
from tempfile import NamedTemporaryFile
import unittest

class TestPil(unittest.TestCase):
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_eps_to_png(self):
        '''Test that converting an EPS file to PNG using PIL works.'''
        src = """%!PS-Adobe 2.0
%%BoundingBox: 0 0 144 144
36 36 72 72 rectfill
/Courier findfont 12 scalefont setfont
36 120 moveto (text) show
showpage
"""
        im = EpsImageFile(StringIO(src))
        im.save(NamedTemporaryFile(suffix='.png'))
