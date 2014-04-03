#!/usr/bin/env python
'''
============================================================
Plot "correlation" matrices of haps of siblings in a family
with non-genotyped parents.   
 
Created on September 7, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import matplotlib.pyplot as P, csv, os, impute as im, networkx as nx

if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    
    p = im.hutt('hutt.npz')
    cm = p.info.snp['dist_cm']
    
    # Region
    a, b = p.info.snp[1050]['dist_cm'], p.info.snp[2150]['dist_cm']

    L = b - a
    G = nx.Graph()
    G.add_weighted_edges_from((i, j, w / L) for (i, j, w) in 
                              (((int(x[4]), int(x[5])), (int(x[6]), int(x[7])),
                                max(0, min(cm[int(x[1])], b) - max(cm[int(x[0])], a))) for x in
                               csv.reader(open(os.environ['OBER_OUT'] + '/ibd/chr22/segments.15cm.out', 'rb'), delimiter=' '))
                              if w > 0)
    c = nx.connected_component_subgraphs(G)
    
    for i in xrange(4):
        P.figure(1)
        P.clf()
        nx.draw(c[i], node_size=100, font_size=6, with_labels=False)
        P.title('# Haplotypes = %d' % (c[i].number_of_nodes(),))
        P.savefig('ibd-graph-%d.png' % (i,))
        P.show()
        
