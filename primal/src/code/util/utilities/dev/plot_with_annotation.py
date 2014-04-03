'''
Created on Dec 6, 2012

@author: oren
'''
import numpy as np
import matplotlib.pyplot as plt

N = 10
data = np.random.random((N, 4))
labels = ['point{0}'.format(i) for i in xrange(N)]

plt.subplots_adjust(bottom=0.1)
plt.scatter(
    data[:, 0], data[:, 1], marker='o', c=data[:, 1], s=100,
    cmap=plt.get_cmap('Spectral'))

for label, x, y in ((label, x, y) for (label, x, y,) in zip(labels, data[:, 0], data[:, 1]) if y < 0.5):
    plt.annotate(
        label,
        xy=(x, y), xytext=(-20, 20),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

plt.show()
