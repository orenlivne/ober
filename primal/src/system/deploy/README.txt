Ober Lab python Code Deployment on Ubuntu (e.g., oberlab-dev)
=============================================================

1. Set up debian repositories
-----------------------------

In  /etc/apt/sources.list:
deb http://debian.uchicago.edu/debian/ squeeze main contrib
deb http://debian.uchicago.edu/debian/ testing main contrib
deb http://debian.uchicago.edu/debian/ unstable main contrib
deb http://security.debian.org/ testing/updates main contrib

sudo apt-get update

2. Install prerequisite debian packages
---------------------------------------

Note: we use C, C++ and GFortran 4.7 versions, which are experimental at this time.
You can probably use a stable version - make sure all three have the same version.

sudo apt-get gcc-4.7 g++-4.7 gfortran-4.7 python python-virtualenv virtualenvwrapper python-mysqldb libmysqlclient-dev python-dev libfreetype6-dev libpng12-0 libpng12-dev subversion

sudo ln -s /usr/bin/gfortran-4.7 /usr/bin/gfortran

3. Install python 2.7 into a directory where you have write-permission.

Here we use /opt/python. The following was done as root since it owns this directory,
but you can do it as your user account under your home directory.

mkdir /opt
wget http://www.python.org/ftp/python/2.7.3/Python-2.7.3.tgz
tar -zxvf Python2.7.3.tar.gz
cd Python2.7.3
ln -s Python2.7.3 python
cd python
./configure --with-zlib --prefix=/opt/python
make
make install


4. Get the Ober Lab Code
------------------------
cd ~
svn https://oberlab-tk.uchicago.edu/svn/ober


5. Set up environment
---------------------

Add the following environment variables to your ~/.bash_profile:

# Directory where you checked out the SVN tree
export OBER="$HOME/ober"
source $OBER/system/dots/bash_profile

source ~/.bash_profile

6. Create new virtual environment "ober"
----------------------------------------

mkdir -p $WORKON_HOME
mkvirtualenv ober -p /opt/python/bin/python

7. Install python packages within the virtual environment
---------------------------------------------------------

workon ober
pip install numpy scipy matplotlib django MySQL-python simplejson networkx nose

Check your versions:
python -c 'import numpy, scipy, matplotlib, django, networkx, simplejson, MySQLdb, nose; print "NumPy version:", numpy.__version__; print "SciPy version:", scipy.__version__; print "MatPlotLib version:", matplotlib.__version__; print "Django version:", django.get_version(); print "NetworkX version", networkx.__version__; print "SimpleJson version:", simplejson.__version__; print "MySQLdb version:", MySQLdb.__version__; print "Nose version:", nose.__version__'

Should print something like that:
NumPy version: 1.6.2
SciPy version: 0.10.1
MatPlotLib version: 1.1.0
Django version: 1.4.1
NetworkX version 1.7
SimpleJson version: 2.6.1
MySQLdb version: 1.2.3


8. Run the hera django webapp
-----------------------------
In this example, we use oberlab-dev's public IP address ( 10.135.144.36;
requires being on the campus VPN to access).

cd $OBER/hera

Drop and recreate database:
./manage.py sqlclear pedtools | ./manage.py dbshell
./manage.py syncdb

./manage.py runserver 10.135.144.36:8000 

Then point your browser to: http://oberlab-dev.uchicago.edu:8000/pedtools/