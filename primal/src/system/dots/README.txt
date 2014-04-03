bash environment scripts for different machines. Include in your ~/.bash_profile to set paths and environment variables useful for running the programs under the "code" SVN directory. Both Linux and Cygwin are supported.

Prerequisite software packages is assumed to be installed under the $APPS directory:
- bash shell (Linux/Cygwin; preferably Linux)
- python tools
  ** python packages: see svn:/ober/system/bin/provision-python for a full list. This script is capable
     of installing all python tools mentioned in this README from scratch.
  ** python 2.7.3 - http://www.python.org/getit/
  ** virtualenv - http://pypi.python.org/pypi/virtualenv
- plink (genetic text data manipulation) - http://pngu.mgh.harvard.edu/~purcell/plink/
- pedfiddler (for pedigree drawing) - http://www.stat.washington.edu/thompson/Genepi/Pedfiddler.shtml
- tabix (for fast exome data lookup) - http://sourceforge.net/projects/samtools/files/tabix/

Development tools:
- virtualenvwrapper - http://pypi.python.org/pypi/virtualenvwrapper
- ipython 0.13 - http://ipython.org/download.html
- Eclipse - http://www.eclipse.org/downloads/
  ** PyDev python development plugin - http://pydev.org/download.html
