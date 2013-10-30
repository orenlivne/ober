#!/usr/bin/env python

def bp_slice(snp, pos, start_bp, end_bp):
     slice = filter(lambda (i,pos):  start_bp < pos < end_bp, enumerate(pos))
     indexes = map(lambda (i,slice):i, slice)
     snpslice = []
     for i in indexes:
          snpslice.append(snp[i])
     print snpslice

     snpslice2 = [snp[i] for i in indexes]
     print snpslice2

     # Numpy indexing would make this even easier!

if __name__ == '__main__':
     snp = ['rs1', 'rs2', 'rs3', 'rs4']
     pos = [100, 200, 300, 400]
     bp_slice(snp, pos, 199, 300)

