'''
============================================================
A (string) comparator example: sort a list of strings
composed of "ACGT?" so that A appears first, then T, then
C, then G, then ?.
============================================================
'''

ORDER = {'a': 0, 't': 1, 'c': 2, 'g': 3, '?': 4} 

def acgt_comparator(s, t):
    s_sig = map(ORDER.__getitem__, s)
    t_sig = map(ORDER.__getitem__, t)
    if s_sig < t_sig: return -1
    elif s_sig == t_sig: return 0
    else: return 1

def acgt_comparator2(s, t):
    k, max_k = 0, min(len(s), len(t))
    while k < max_k and ORDER[s[k]] == ORDER[t[k]]: k += 1
    if k == max_k: return -1 if len(s) < len(t) else 1  # Both substrings equal to end of string so shorter one is smaller
    else: return -1 if ORDER[s[k]] < ORDER[t[k]] else 1  # Unequal at a character before end of string
    
if __name__ == "__main__":
    a = ['acgt', 'cagt', 'tggg', 'cag', 'gat', 'ga?', 'ga']
    print sorted(a, cmp=acgt_comparator)
    print sorted(a, cmp=acgt_comparator2)
