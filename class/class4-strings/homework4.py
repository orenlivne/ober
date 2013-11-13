'''
============================================================
Programming Pearls, Section 15.1: counting words.
============================================================
'''
from __future__ import division
import itertools as it, re
from math import ceil
from collections import Counter
import cProfile, pstats

WORD_PATTERN = re.compile('^([a-zA-Z]*)((\-[a-zA-Z]*)*)[\?|\!]*$')

def to_word(s):
    m = WORD_PATTERN.match(s.lower())
    if m:
        groups = m.groups()
        yield groups[0]
        if len(groups) >= 1:
            for w in groups[1][1:].split('-'): yield w

words = lambda f: it.ifilter(lambda x: x, (word for line in f for s in line.split() for word in to_word(s)))

def column_format(entries, cols, delimiter='\t'):
    n = len(entries)
    col_size = int(ceil(n / cols))
    for j in xrange(col_size):
        for i in xrange(j, n, col_size):
            if i > j: print '%s' % (delimiter,),
            w, c = entries[i]
            print '%-10s %-5d' % (w, c),
        print ''
    
def print_top_entries(file_name, entries=21, cols=3):
    with open(file_name, 'rb') as f: column_format(Counter(words(f)).most_common(entries), 3)

def print_all_entries(file_name, cols=3):
    with open(file_name, 'rb') as f:
        for w, c in Counter(words(f)).iteritems(): print w, c

if __name__ == "__main__":
    cProfile.run('print_top_entries("bible.txt", entries=200)', 'homework4')
    p = pstats.Stats('homework4').strip_dirs()
    p.sort_stats('cumulative').print_stats(50)
