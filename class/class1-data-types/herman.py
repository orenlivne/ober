snps = ['rs2980300','rs6685064','rs3766180','rs6603791','rs7519837','rs2281173','rs1107910','rs2272908','rs3737628','rs9786963','rs10907187','rs7511905']
bps = [785989,1211292,1478153,1500941,1510801,1688192,1692321,1721479,1722932,1759026,1759054,1793786]


def get_snps(snp,bp,start,stop):
    selected = {}
    for rs, position in zip(snps, bps):
        if position >= start and position <= stop: 
            selected[rs] = position
    return selected

print get_snps(snps,bps,70000,1200000)
