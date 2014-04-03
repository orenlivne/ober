'''
Created on Jan 17, 2013

Print original graph weights and and affinities for edges between sets of conflicting alleles (R1,R2)
identified at chr22:SNP#456 imputation.

@author: oren
'''
import pickle, os

def load_graph(snp):
    print 'load_graph(%d)' % (snp,)
    return pickle.load(open(os.environ['OBER_OUT'] + '/ibd/chr22/index_segments_%d/graph-%d.pickle' % (snp, snp), 'rb'))

G = load_graph(457)
