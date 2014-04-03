# Subjects 80, 248 are sibs in the pedigree but might be only half-sibs. Same mom different dads.
# Try to confirm/reject this hypothesis and locate the other dad based on IBD segments.

import numpy as np, util, impute as im, matplotlib.pyplot as P
from impute.data.constants import PATERNAL

def autolabel(ax, rects):
    '''Attach some text labels.'''
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, '%d' % int(height),
                ha='center', va='bottom')

def closest_relatives(sample, num=50):
    '''IBD sharing information on num closest relatives.'''
    s = np.loadtxt('/home/oren/ober/out/validation/affy/segments.%d.%d.out' % (sample, PATERNAL), dtype=int)
    values, groups = util.sum_by_group(s[:, 3] - s[:, 2], 2 * s[:, 4] + s[:, 5])  # hash haps columns 4 and 5
    j = np.argsort(values)[-1::-1][:num]
    return np.concatenate((values[j][np.newaxis].transpose() / total,
                           groups[j][np.newaxis].transpose(),
                           np.array([im.pt.lowest_common_ancestor(ped.graph, sample, relative) for relative in groups[j] / 2])),
                          axis=1)
    
def bar_data(r0, r1, num=50):
    '''Convert relative data to two-group bar chart data.'''
    def to_key(r):
        x = r[:, 1].astype(int)
        return [(y / 2, y % 2) for y in x]
    d0 = dict(zip(to_key(r0), [(x[0], int(x[3])) for x in r0]))
    d1 = dict(zip(to_key(r1), [(x[0], int(x[3])) for x in r1]))

    d = {}
    for k in set(d0.keys()) | set(d1.keys()):
        d.setdefault(k, [(0, 0), (0, 0)])
        if k in d0: d[k][0] = d0[k]
        if k in d1: d[k][1] = d1[k]
    for k in d:
        v0, v1 = d[k]
        d[k] = [v0[0], v1[0], max(v0[1], v0[1])]
        
    entries = sorted(d.iteritems(), key=lambda (k, v):-v[1])[:num]
    for x in entries:
        print x
    return [x[1][0] for x in entries], \
        [x[1][1] for x in entries], [(x[0], x[1][2]) for x in entries], d

def plot_ibd_sharing(k0, k1, labels):
    N = len(k0)

    ind = np.arange(N)  # the x locations for the groups
    width = 0.35  # the width of the bars
    
    P.close(1)
    P.figure(num=1, figsize=(15, 6))
    P.clf()
    ax = P.axes()
    rects1 = ax.bar(ind, k0, width=width, color='b')
    rects2 = ax.bar(ind + width, k1, width=width, color='r')
    
    # add some
    ax.set_ylabel('Genomic Fraction of IBD Sharing')
    ax.set_title('Sample IBD Sharing Comparison')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(['(%d,%d)\n%d' % (x[0][0], x[0][1], x[1]) for x in labels])
    ax.legend((rects1[0], rects2[0]), ('Sample 80', 'Sample 248'))
#    autolabel(ax, rects1)
#    autolabel(ax, rects2)
    P.show()

####################################################################################
if __name__ == '__main__':
    total = sum(np.loadtxt('/home/oren/ober/testdata/misc/chromosome.txt', usecols=[3]))
    ped = im.hutt_pedigree()
    np.set_printoptions(linewidth=100)
    samples = (80, 248)
    r0 = closest_relatives(80)
    r1 = closest_relatives(248)
    
    k0, k1, labels, d = bar_data(r0, r1, num=50)
    plot_ibd_sharing(k0, k1, labels)
    P.savefig('/home/oren/ober/doc/imputation/validation/affy/ibd_comparison.png')
    
    k0, k1, labels, d = bar_data(r0, r1, num=8)
    plot_ibd_sharing(k0, k1, labels)
    P.savefig('/home/oren/ober/doc/imputation/validation/affy/ibd_comparison_top.png')
        