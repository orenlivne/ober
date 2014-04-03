#!/usr/bin/env python
'''
============================================================
Manages an SQLITE3 SNP annotation database. 

Created on October 17, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, csv, networkx as nx, numpy as np, util
from optparse import OptionParser
from db_gene.snp.entities import Chromosome, Base, Snp, Ld, Kinship
from db_gene.ucsc import ucsc_dao
from db_gene import snp
from sqlalchemy import func
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from db_gene.snp.snp_db_dao import Daos

#---------------------------------------------
# Private Methods
#---------------------------------------------
def __read_chrom_metadata(file_name):
    '''Read chromosome length metadata from the file file_name into a record array.'''
    return np.loadtxt(file_name,
                      usecols=range(0, 5),
                      dtype={'names': ('number',
                                       'name',
                                       'num_genes',
                                       'total_bp',
                                       'sequenced_bp'),
                             'formats': (np.int, np.chararray, np.int, np.long, np.long)})

def __read_ld_edges(in_file, threshold, report_once_in=100000):
    '''Read LD edges (i,j,LD[i,j]) from the input file in_file. Filter entries smaller than threshold.'''
    num_records = 0
    for (k, line) in enumerate(csv.reader(in_file, delimiter=' ', skipinitialspace=True)):
        if num_records and np.mod(k, report_once_in) == 0:
            print 'Read %d lines, %d #records so far' % (k, num_records)
        if line:
            (i, j, w) = (line[3], line[4], float(line[6]))
            if w >= threshold:
                num_records += 1
                yield (i, j, w)

def __read_map(in_file, report_once_in=100000, max_records=None):
    '''Read a chromosome SNP genetic map file. Yield (chrom_number, SNP name, base pair position).'''
    num_header_lines = 1
    reader = csv.reader(in_file, delimiter='\t', skipinitialspace=True)
    for line in xrange(0, num_header_lines): reader.next()
    num_records = 0
    for (k, line) in enumerate(reader):
        if num_records and np.mod(k, report_once_in) == 0:
            print 'Read %d lines, %d #records so far' % (k, num_records)
        if line:
            num_records += 1
            yield (int(line[0]), line[1], line[2])
        if num_records == max_records:
            break

def read_snp_data(bim_file, snp_dao):
    '''Read Hutterites SNP data from a PLINK BIM (binary map) file and save to our SNP table.'''
    sys.stdout.write('Loading BIM data ...\n')
    sys.stdout.flush()  
    data = np.loadtxt(bim_file, usecols=[0, 1, 3], dtype=[('chrom', 'i4'), ('name', 'S20'), ('bp', 'i8')])
    sys.stdout.write('Converting to entities ...\n')
    sys.stdout.flush()  
    entities = [snp.entities.Snp(row['chrom'], row['name'], row['bp']) for row in data]
    sys.stdout.write('Persisting %d entities ...\n' % (len(entities),))
    sys.stdout.flush()  
    snp_dao.save(entities)
    sys.stdout.write('Total SNP records: %d\n' % snp_dao.num_records())
    sys.stdout.flush()  

def calculate_genetic_positions(map_file, snps):
    '''Calculate the closest genetic positions to a list of SNP bp''s on a chromosome in the
    chromosomal genetic map file map_file, and update the SNP database with those positions.'''
    positions = np.loadtxt(map_file, skiprows=1, usecols=[1, 3], dtype=np.float)
    bp, genetic = positions[:, 0], positions[:, 1]
    closest_bp = list(util.nearest_neighbor_in_sorted_arrays(bp, [x.bp for x in snps]))
    for i, x in enumerate(snps): x.genetic_pos = genetic[closest_bp[i]]
    snp_dao.save(snps)

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse and validate command-line arguments
    PROGRAM = os.path.basename(sys.argv[0])
    usage = 'Usage: %s <db-url> <command> [command-args]\n\n' \
        'Manage a genomic annotation database.\n\n' \
        'Commands:\n' \
        '\tcreate - Create a new empty database schema\n' \
        '\treport - Print a report with # items of each type\n' \
        '\timport_chr - Import chromosome metadata from text file into the db\n' \
        '\timport_ld - Download LD information from a HapMap file name/URL import it into the db\n' \
        '\timport_snps - Read SNPs from a NOF file, fetch their info from the UCSC table and save in the SNP database table\n' \
        '\timport_genetic_dist - Read SNP genetic positions from genetic map file to SNP db table\n' \
        '\nType ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = OptionParser(usage=usage)
    parser.add_option('-v', '--debug'          , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
#    parser.add_option('-s', '--stage'          , type='int'           , dest='stage', default=0,
#                      help='Run only a this phasing stage')
#    parser.add_option('-f', '--impute', type='int', dest='impute', default=phase.IMPUTE_OPTION.NONE,
#                      help='Post-processing: do nothing (0), impute genotypes from called haplotypes \
#                      (1), or impute and fill missing genotypes randomly from estimated frequencies (2)\
#                      (default)')
#    parser.add_option('-z', '--zero-partial'   , action='store_true'  , dest='zero_partial', default=True,
#                      help='Post-processing: zero-out partially-called genotypes')
    (options, args) = parser.parse_args(sys.argv[1:])
    options.print_times = True
    if len(args) < 2:
        print usage
        sys.exit(1)    
    db_url = args[0]
    command = args[1]

    #db_url = 'sqlite:///%s' % (db_url,)
    engine = create_engine(db_url, echo=options.debug)
    Session = sessionmaker(bind=engine)
    session = Session()
    ld_dao = Daos(db_url, echo=options.debug).ld_dao
    snp_dao = Daos(db_url, echo=options.debug).snp_dao
    ucsc_snp_dao = ucsc_dao.Daos(db_url, echo=options.debug).snp_dao
    
    if command == 'create':
        # Create a new database schema in db_url
        print 'Creating database'
        Base.metadata.create_all(engine)
    elif command == 'report':
        # Show # items of each type
        print 'Database Report'
        print '---------------'
        print 'Chromosome\t: %d' % (session.query(func.count(Chromosome.id)).scalar(),)
        print 'SNP\t\t: %d' % (session.query(func.count(Snp.id)).scalar(),)
        print 'LD\t\t: %d' % (session.query(func.count(Ld.id)).scalar(),)
        print 'Kinship\t\t: %d' % (session.query(func.count(Kinship.id)).scalar(),)
    elif command == 'import_chr':
        # Import chromosome metadata from text file into the database
        if len(args) != 3:
            print usage
            print '\nSyntax: %s <file>' % (command,)
            sys.exit(1)
        # Load data as a record array
        records = __read_chrom_metadata(args[2])
        # Delete existing data
        session.query(Chromosome).delete()
        session.commit()
        # Convert to an entity list and persist to db
        session.add_all(Chromosome(r['number'], r['name'], r['num_genes'], r['total_bp'], r['sequenced_bp']) for r in records)
        session.commit()
        print 'Imported %d chromosome records.' % (len(records),)
#    elif command == 'import_map': # Very slow; better to directly load data from UCSC into sqlite3 if at all needed
#        # Import the genetic map of a single chromosome
#        if len(args) != 3:
#            print usage
#            print '\nSyntax: %s <URL>' % (command,)
#            sys.exit(1)
#        url = args[2]
#        # Read data from a genetic map generated by a UCSC database query to the snp<build_num> table,
#        # convert to entities and persist to db
#        session.add_all(Snp(chrom, name, bp) for (chrom, name, bp) in __read_map(open_resource(url)))
#        session.commit()
    elif command == 'import_ld':
        # Download latest LD information from HapMap website and save it in the database.
        # Compute SNP module numbers and save them as well.
        if len(args) != 4:
            print usage
            print '\nSyntax: %s <chrom> <threshold>' % (command,)
            sys.exit(1)
        chrom = int(args[2])
        threshold = float(args[3])

        # Read LD data
        print 'Loading and building LD graph from db ...'
        g = ld_dao.ld_graph(chrom, threshold=threshold)
        
        # Calculate SNP module indices = connected component index of the corresponding g-nodes.
        # Indices are 0-based.  
        print 'Calculating modules ...'
        n = g.number_of_nodes()
        module = np.zeros((n,), dtype=np.uint)
        for (i, component) in enumerate(nx.connected_components(g)):
            module[component] = i
    elif command == 'import_snps':
        # Read SNPs from a NOF file, fetch their info from the UCSC table and save in the SNP
        # database table
        if len(args) != 3:
            print usage
            print '\nSyntax: %s <bim_file>' % (command,)
            sys.exit(1)
        bim_file = args[2]

        print 'Loading SNPs ...'
        g = read_snp_data(bim_file, snp_dao)
    elif command == 'import_genetic_dist':
        # Load genetic distances from a genetic map file into SNP table
        if len(args) != 4:
            print usage
            print '\nSyntax: %s <chrom> <map_file>' % (command,)
            sys.exit(1)
        chrom = int(args[2])
        map_file = args[3]

        snps = snp_dao.get(chrom)
        calculate_genetic_positions(map_file, snps)
    else:
        print usage
        print '\nUnrecognized command: ''%s''' % (command,)
        sys.exit(1)

    # Clean up        
    session.close()
