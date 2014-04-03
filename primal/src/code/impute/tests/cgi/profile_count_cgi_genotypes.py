#!/usr/bin/env python
# encoding: utf-8
# filename: profile_count_cgi_genotypes.py

import pyximport  # @UnresolvedImport
pyximport.install()

from impute.cgi import count_cgi_genotypes
import pstats, cProfile, os

def count():
    count_cgi_genotypes.main(data_file=os.environ['OBER_OUT'] + '/impute_cgi/imputed_cgi.chr22.head.tsv',
                             variant_id=True)

cProfile.run('count()', 'count')
s = pstats.Stats('count')
s.strip_dirs().sort_stats('time').print_stats()
