'''
============================================================
HUTT Genome browser - gene annotation Data Access Object (DAO).

Created on November 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import numpy as np
from db_gene.dao import Dao  # Must be the first import
from sqlalchemy import func
from sqlalchemy.sql.expression import and_
from db_gene import DEFAULT_URL
from annotate.entities import Variant
from collections import Counter
    
'''Standardized category names in standardized order.'''
REGION_CATEGORIES = ['Coding', 'Splicing', 'Intronic', '3UTR', '5UTR + TSS Upstream', 'Intergenic']

'''Maps variant category annotation column to standardized category name.''' 
REGION_CATEGORY_OF = dict([
      ('-', 'Intergenic'),
      ('downstream', '5''UTR + TSS Upstream'),
      ('exonic', 'Coding'),
      ('exonic;splicing', 'Splicing'),
      ('intergenic', 'Intergenic'),
      ('intronic', 'Intronic'),
      ('upstream', '5''UTR + TSS Upstream'),
      ('ncRNA_exonic', 'Coding'),
      ('ncRNA_intronic', 'Intronic'),
      ('ncRNA_splicing', 'Splicing'),
      ('ncRNA_UTR3', '3''UTR'),
      ('ncRNA_UTR5', '5''UTR + TSS Upstream'),
      ('splicing', 'Splicing'),
      ('upstream;downstream', '5''UTR + TSS Upstream'),
      ('UTR3', '3''UTR'),
      ('UTR5', '5''UTR + TSS Upstream'),
      ('UTR5;UTR3', '5''UTR + TSS Upstream')
      ])
    
'''Standardized coding variant category names in standardized order.'''
CODING_CATEGORIES = ['Missense+start+stop', 'Nonsense', 'Synonymous', 'Frameshift', 'In-frame Indels']

'''Maps variant coding variant category annotation column to standardized category name.''' 
CODING_CATEGORY_OF = dict([
      ('-', None),
      ('DELETE', 'In-frame Indels'),
      ('DELETE+', 'Frameshift'),
      ('DISRUPT', None),
      ('FRAMESHIFT', 'Frameshift'),
      ('INSERT', 'In-frame Indels'),
      ('INSERT+', 'Frameshift'),
      ('MISSENSE', 'Missense + Misstart + Misstop'),
      ('MISSTART', 'Missense + Misstart + Misstop'),
      ('NO-CHANGE', 'Synonymous'),
      ('NONSENSE', 'Nonsense'),
      ('NONSTOP', 'Missense + Misstart + Misstop'),
      ('SYNONYMOUS', 'Synonymous'),
      ('UNKNOWN-INC', None),
      ('UNKNOWN-TR', None)
      ])

####################################################################################
class VariantSummaryReport(object):
    '''Represents a variant count summary report, broken down by variant category and MAF (binned).'''
    def __init__(self, categories, maf, count):
        '''Columns correspond to the label list ''categories'', column i to the MAF bins
        (maf[i],maf[i+1]].''' 
        self._categories, self._maf, self._count = categories, maf, count

    @staticmethod
    def load(input_stream, delimiter='\t'):
        '''Load report data from input stream ''input''.'''
        categories = next(input_stream).strip().split(delimiter)[1:]
        data = np.array([map(float, line.strip().split(delimiter)) for line in input_stream])
        maf, count = np.concatenate(([data[0, 0]], data[:, 1])), data[:, 2:].astype(int)
        return VariantSummaryReport(categories, maf, count)
        
    def save(self, out, delimiter='\t'):
        '''Save report to the output stream ''out''.'''
        out.write('MAF' + delimiter + delimiter.join('%s' % (category,) for category in self._categories) + '\n')
        for i in xrange(len(self._maf) - 1):
            out.write('%f%s%f' % (self._maf[i], delimiter, self._maf[i + 1]) + delimiter + delimiter.join('%d' % (x,) for x in self._count[i, :]) + '\n')

####################################################################################
class VariantDao(Dao):
    '''Data access object for variant annotations.'''
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        super(VariantDao, self).__init__(url, **kwargs)
    
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def num_records(self):
        '''Return the total number of records in the refVariant database.'''
        # return self.engine.execute('select count(*) from refVariant').scalar()
        session = self.Session()
        result = self.Session().query(func.count(Variant.record_id)).scalar()
        session.close()
        return result

    def variant_count_report_maf_bin(self, group, maf_min, maf_max):
        '''Return variant count for a MAF range (maf_min, maf_max], broken down by region functional annotations.'''
        # print 'intersecting_genes', chrom, start, end
        _, _, field = VariantDao.group_abstract_factory(group) # Appear in standardized order in result object
        session = self.Session()
        result = list(session.query(field, func.count()).\
                      filter(and_(Variant.maf_imputed > '%f' % (maf_min,),
                                  Variant.maf_imputed <= '%f' % (maf_max,),
                                  Variant.is_qc)).\
                      group_by(field))
#        print result
        session.close()
        return result

    def variant_count_report(self, group, maf=None):
        '''Generate a variant region count report. Uses default MAF bins unless a non-null numpy array
        is specified in the ''maf'' property.'''
        categories, category_of, _ = VariantDao.group_abstract_factory(group) # Appear in standardized order in result object
        maf = maf if maf is not None else np.concatenate((np.arange(0, 0.11, 0.01), [0.5]))
        count = np.zeros((len(maf) - 1, len(categories)), dtype=int)
        for i in xrange(len(maf) - 1):
            result = Counter()
            for k, v in ((category_of[k], v) for k, v in self.variant_count_report_maf_bin(group, maf[i], maf[i + 1])): result[k] += v
            count[i, :] = [result[category] for category in categories]
        return VariantSummaryReport(categories, maf, count)
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------
    @staticmethod
    def group_abstract_factory(group):
        if group == 'region':
            return REGION_CATEGORIES, REGION_CATEGORY_OF, Variant.var_region
        elif group == 'coding':
            return CODING_CATEGORIES, CODING_CATEGORY_OF, Variant.var_mutation
        else:
            raise ValueError('Unrecognized group ''%s''' % (group,))

####################################################################################
class Daos(Dao):
    '''DAO mother object.'''
    def __init__(self, url, **kwargs):
        self.configure(url, **kwargs)

    def configure(self, url, **kwargs):
        # Synchronize access to url
        self.variant_dao = VariantDao(url, **kwargs)

# Global access to DAO, default configuration. Lazily-initialized, to accommodate systems without mysql
__DEFAULT_HUTT_DAOS = None

def DEFAULT_HUTT_DAOS():
    global __DEFAULT_HUTT_DAOS
    if not __DEFAULT_HUTT_DAOS:
        __DEFAULT_HUTT_DAOS = Daos(DEFAULT_URL)
    return __DEFAULT_HUTT_DAOS
