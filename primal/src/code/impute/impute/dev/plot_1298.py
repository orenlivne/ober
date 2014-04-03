#!/usr/bin/env python
'''
============================================================
Plot the nbhrs1298 pedigree.

Created on August 16, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, famplot as fp

p = im.io.read_npz(im.itu.NBHRS1298_STAGE4)
in_file = 'NBHRS1298_STAGE4.txt'
im.pt.to_pedfiddler(p.pedigree, in_file)

coord_params = fp.CoordParams()
coord_params.algorithm = 'default'

out = 'NBHRS1298_STAGE4.eps'
# out = tempfile.NamedTemporaryFile(suffix='.eps', delete=False)
ped_info = fp.draw_pedigree(in_file, coord_params,
                        fp.PlotParams(), fp.DEFAULT_COLORS, out)
