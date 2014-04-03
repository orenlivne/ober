#!/usr/bin/env python
'''
============================================================
Plot the nbhrs1298 pedigree.

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im

p = im.hutt('hutt.phased.npz')

s = 1159  # Unphased sample index
#phaser = im.phase_core.new_phaser_chain([im.phase_distant.distant_phaser(single_sample=s)])
phaser = im.phase_core.new_phaser_chain([im.phase_distant.distant_phaser(phased_fill=0.92, target_fill=0.95, max_path_length=7)])
h = p.haplotype
print h.fill_fraction(sample=s)
phaser.run(p, im.PhaseParam(debug=True))
print h.fill_fraction(sample=s)

out = '/home/oren/ped-%d.png' % (s,)
im.pt.draw_member_neighbor_genotyped_pedigree(p, 1159, 5, out, identifier='index')
