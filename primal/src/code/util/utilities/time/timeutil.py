'''
============================================================
Time utilities.

Created on April 3, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import datetime, time

def utc2datetime(s, fmt='%a %b %d %H:%M:%S %Y'):
    '''Convert a string like ''Wed Apr  3 6:18:03 2013'' to a datetime object.''' 
    return datetime.datetime.fromtimestamp(time.mktime(time.strptime(s, fmt)))

def str2timedelta(s):
    '''Convert a string in the format [DD:]HH:MM:SS to a datetime.timedelta object.'''
    parts = map(int, s.split(':'))
    if len(parts) > 0:
        # Pad days if not present
        parts.insert(0, 0)
    return datetime.timedelta(days=parts[0], hours=parts[1], minutes=parts[2], seconds=parts[3])

def format_timedelta(delay, pad_days=False):
    out = str(delay)
    if delay.days > 0:
        out = out.replace(" days, ", ":")
    elif pad_days:
        out = "0:" + str(delay)
    outAr = out.split(':')
    outAr = ["%02d" % (int(float(x))) for x in outAr]
    out = ":".join(outAr)
    return out
