import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext  # @UnresolvedImport

ext_modules = [Extension("impute.cgi.count_cgi_genotypes2", ["impute/cgi/count_cgi_genotypes2.pyx"])]

# Utility function to read the README file.  
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='impute',
    version='0.0.1',
    author='Oren Livne',
    author_email='livne@uchicago.edu',
    description=('Gene imputation tools for founder populations.'),
    license='BSD',
    keywords='gene impute',
    url='http://packages.python.org/impute',
    packages=['impute', 'tests'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    long_description=read('README'),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
        'License :: OSI Approved :: BSD License',
    ],
)
