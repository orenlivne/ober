#!/usr/bin/env python
'''
============================================================
Convert coordinates read from stdin from one human genome
build to another.

Usage: lift_over.py <from-build> <to-build>

stdin line format: chrom bp_in_from_build
stdout line format: bp_in_to_build, or '-' if not found

Created on February 19, 2014
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, traceback, util
from pyliftover import LiftOver

if __name__ == '__main__':
    try:
        src, target = sys.argv[1:3]
        if src == target:
        	for _, bp in (line.strip().split(' ') for line in sys.stdin):
        	    print '%d %d' % (int(bp), int(bp))
        else:
            lo = LiftOver(src, target)
            for chrom, bp in (line.strip().split(' ') for line in sys.stdin):
                out = lo.convert_coordinate('chr' + chrom, int(bp))
                if not out:
                    print '-'
                else:
                    print '%d' % (out[0][1],)
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
