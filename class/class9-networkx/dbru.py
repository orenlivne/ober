'''
============================================================
http://rosalind.info/problems/dbru

Given: A collection of up to 1000 DNA strings of equal length (not exceeding 50 bp) corresponding to a set S of (k+1)-mers.

Return: The adjacency list corresponding to the de Bruijn graph corresponding to S U Src.
============================================================
'''
import itertools as it, networkx as nx

read_lines = lambda file_name: [x.strip() for x in open(file_name, 'rb')]
COMPLEMENT = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
revc = lambda s: ''.join(map(COMPLEMENT.get, reversed(s)))

de_bruijn_adj_list = lambda S: set((r[:-1], r[1:]) for r in it.chain(S, it.imap(revc, S)))
de_bruijn_graph = lambda r: nx.from_edgelist(de_bruijn_adj_list(r), create_using=nx.DiGraph())

if __name__ == "__main__":
#    print dbru('rosalind_dbru_sample.dat')   
    reads = read_lines('/home/oren/class/ober/class/class9-networkx/reads.txt')
    G = de_bruijn_graph(reads)
    nx.draw(G)
    
