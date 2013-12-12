'''
============================================================
http://rosalind.info/problems/lcsq

Given: Two DNA strings s and t (each having length at most 1 kbp) in FASTA format.

Return: A longest common subsequence of s and t. (If more than one solution exists, you may return any one.)
============================================================
'''
import numpy as np

def lcs_length_table(x, y):
    '''DP to find LCS length of all sub-problems.'''
    m, n = len(x), len(y)
    c = np.zeros((m + 1, n + 1), dtype=int)
    for i in xrange(1, m + 1):
        for j in xrange(1, n + 1): c[i, j] = (c[i - 1, j - 1] + 1) if x[i - 1] == y[j - 1] else max(c[i - 1, j], c[i, j - 1])
    return c

def lcsq(x, y):
    '''Back-track to find LCS.'''
    m, n, c = len(x), len(y), lcs_length_table(x, y)
    i, j, s = m, n, ''
    while i > 0 and j > 0:
        if x[i - 1] == y[j - 1]:  i, j, s = i - 1, j - 1, x[i - 1] + s
        else:
            if c[i - 1, j] > c[i, j - 1]: i -= 1
            else: j -= 1
    return s

def lcsq2(x, y):  # Integrated DP+backtracking, O(mn) time, O(min(m,n)) storage
    m, n = len(x), len(y)
    if m < n: return lcsq2(y, x)
    c_old, c, s_old, s = [0] * (n + 1), [0] * (n + 1), [''] * (n + 1), [''] * (n + 1)
    for xi in x:
        c_old[:] = c[:]; s_old[:] = s[:]
        for j, yj in enumerate(y, 1):
            if xi == yj: c[j], s[j] = c_old[j - 1] + 1, s_old[j - 1] + xi
            else:
                if c_old[j] > c[j - 1]: c[j], s[j] = c_old[j], s_old[j]
                else: c[j], s[j] = c[j - 1], s[j - 1]
    return s[-1]

if __name__ == "__main__":
    print lcsq('TAC', 'TCA'), lcsq2('TAC', 'TCA')
    print lcsq('AACCTTGG', 'ACACTGTGA'), lcsq2('AACCTTGG', 'ACACTGTGA')