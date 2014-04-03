import os, sys, site

# Add the site-packages of the chosen virtualenv to work with
VENV='/home/oren/virtualenvs/ober'
site.addsitedir(VENV + '/lib/python2.7/site-packages')

# Add the app's directory to the PYTHONPATH
sys.path.insert(0, '/home/oren/ober/code/annotate')

# Activate your virtual env
activate_env=os.path.expanduser(VENV + '/bin/activate_this.py')
execfile(activate_env, dict(__file__=activate_env))

from flaskr import app as application
