'''
============================================================
Hutterities WGS Data (CGI) ID mappings.

Created on November 28, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import csv, os, numpy as np
from collections import OrderedDict

#---------------------------------------------
# Constants
#---------------------------------------------
'''Location of CGI-to-sample-id index file.'''        
DEFAULT_ID_INDEX = os.environ['OBER_DATA'] + '/cgi/README.assembly_sample_subject.csv'

#---------------------------------------------
# Methods
#---------------------------------------------
def cgi_id_to_sample_id(file_name=DEFAULT_ID_INDEX):
    '''Return a dictionary of WGS Hutterities CGI-ID to FINDIV (Sample ID) read from the CGI
    README summary file file_name.'''             
    return OrderedDict((line[0], int(line[2])) for i, line in
                       enumerate(csv.reader(open(file_name, 'rb'), delimiter=',', skipinitialspace=True))
                       if (i >= 1) and line)

def sample_ids(file_name=DEFAULT_ID_INDEX, sort_ids=False):
    '''Return the list of WGS Hutterities FINDIVs (sample ids).'''
    return np.array(sorted(cgi_id_to_sample_id(file_name).values())) if sort_ids else np.array(cgi_id_to_sample_id(file_name).values())
