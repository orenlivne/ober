'''
Created on Jan 17, 2013

Test-vector edges

@author: oren
'''
import impute as im, pickle, networkx as nx, matplotlib.pyplot as P, os
# from scipy.sparse.csr import csr_matrix

G = pickle.load(open(im.itu.SEGMENT_INDEX_CHR22 + '/graph-1200.txt', 'rb'))
G = nx.connected_component_subgraphs(G)[0]

out_dir = os.environ['OBER'] + '/doc/ibd/index-segments'
P.figure(1)
P.clf()
nx.draw(G)
P.savefig(out_dir + '/G.png')

c, G2, cliques = im.index.partition_amg.partition(G, theta=0.995, full_output=True)
print len(nx.connected_components(G2))

P.figure(2)
P.clf()
P.hist(c, bins=100, log=True)
P.savefig(out_dir + '/affinities.png')

P.figure(3)
P.clf()
nx.draw(G2)
P.savefig(out_dir + '/G-trimmed.png')
