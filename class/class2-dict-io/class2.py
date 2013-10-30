#!/usr/bin/env python
# Example of filtering a SNP list by a filter on a corresponding base-pair list.

def in_range(snp, pos, start_bp, end_bp):
    filtered = []
    for i in xrange(len(snp)):
        if start_bp <= pos[i] <= end_bp:
            filtered.append(snp[i])
    return filtered

def in_range2(snp, pos, start_bp, end_bp):
    filtered = []
    for cur_snp, cur_pos in zip(snp, pos):
        if start_bp <= cur_pos and cur_pos <= end_bp:
            filtered.append(cur_snp)
    return filtered

def in_range3(snp, pos, start_bp, end_bp):
    return [cur_snp for (cur_snp, cur_pos) in zip(snp, pos) 
                if start_bp <= cur_pos and cur_pos <= end_bp]

def in_range4(snp, pos, start_bp, end_bp):
    return (cur_snp for (cur_snp, cur_pos) in zip(snp, pos) 
                if start_bp <= cur_pos and cur_pos <= end_bp)

def read_file(file_name):
    '''Read a snp name column and a bp column from the file name f.'''
    f = open(file_name, 'rb')
    snp = []
    bp = []
    for line in (line.strip() for line in f):
        items = line.split()       
#         '12345 string' --> ['12345', 'string']
        snp.append(items[0])
        bp.append(int(items[1]))
    
    f.close()
    return snp, bp

def read_file2(file_name):
    '''Read a snp name column and a bp column from the file name f.'''
    with open(file_name, 'rb') as f:
        return zip(*(((items[0], int(items[1])) for items in (line.strip().split() for line in f))))    

if __name__ == '__main__':
    snp = ['rs2980300', 'rs6685064', 'rs3766180', 'rs6603791', 'rs7519837', 'rs2281173', 'rs1107910', 'rs2272908', 'rs3737628', 'rs9786963', 'rs10907187', 'rs7511905']
    pos = [785989, 1211292, 1478153, 1500941, 1510801, 1688192, 1692321, 1721479, 1722932, 1759026, 1759054, 1793786]
    print in_range(snp, pos, 70000, 1400000)
    print in_range2(snp, pos, 70000, 1400000)
    print in_range3(snp, pos, 70000, 1400000)
    print in_range4(snp, pos, 70000, 1400000)
    
    filtered = in_range3(snp, pos, 70000, 1400000)
    for cur_snp in filtered:
        print '%-10s' % (cur_snp)

    for cur_snp in filtered:
        print '%-10s' % (cur_snp)

    x, bp = read_file2('/home/oren/ober/code/misc/class/input.txt')
    print 'snp', x
    print 'bp', bp
