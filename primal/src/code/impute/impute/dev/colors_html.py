'''
============================================================
Show our color palette in an HTML document.

Created on Jan 14, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, os

out_file = os.environ['OBER'] + '/doc/ibd/colors.html'

with open(out_file, 'wb') as out:
    out.write('<html><body>\n')
    out.write('\n'.join(map(lambda x: '<div style="width:30px; height: 10px; background-color:%s;float: left"> </div>' % (x,), \
                        im.plot.colors.gethtmlcolors(100))))
    out.write('</body></html>\n')
