#!/usr/bin/env python
'''
============================================================
Cluster job time display (start, stop, projected end).
Input is read from stdin:

job start time (ISO format)
job walltime ([DD:]HH:MM:SS)
job completion fraction (between 0 and 1)

Output:
remaining time
projected remaining time

Created on April 3, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import traceback, sys, util, datetime
from utilities.time import timeutil

#---------------------------------------------
# Methods
#---------------------------------------------
'''Round a timedelta object to the nearest second.'''
round_to_seconds = lambda x: datetime.timedelta(seconds=round(x))

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    try:
        args = [x.rstrip('\n').strip() for x in sys.stdin.readlines()]
        #print args
        start = timeutil.utc2datetime(args[0])
        now = datetime.datetime.now()
        wallclock = timeutil.str2timedelta(args[1])
        f = max(0.001,float(args[2]))  # Job completion fraction (in [0,1])

        remaining = round_to_seconds((start + wallclock - now).total_seconds())
        projected_remaining = round_to_seconds((now-start).total_seconds()*((1.-f)/f))

        print str(remaining), str(projected_remaining), int(remaining > projected_remaining)
        #print projected.strftime('%a %b %d %H:%M:%S %Y')
        
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
