'''
============================================================
Simulate gene dropping (only approximately, via the
distribution of detailed identity states) to estimate the
upper bound on the full- and partial-genotype call rates.

Version 1 - doesn't work well because it ignores the IBD
structure among the 98, leading to inflated call rates.
Instead, use gene dropping simulations (call_rate_sim.py).

Created on March 7, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import impute as im, numpy as np, statutil, os, sys, optparse, util, traceback

#---------------------------------------------
# Constants
#---------------------------------------------

# Recoding detailed identity states to A-allele imputability (A0=left bit; A1=right bit)
IMPUTABLE = np.array([3, 3, 3, 2, 1, 0, 0, 0, 3, 2, 1, 3, 2, 1, 0])

# Probability of a detailed sub-state within its condensed state
SUB_PROB = np.array([1, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 0.5, 0.25, 0.25, 0.5, 0.25, 0.25, 1])

# Condensed state identifier of detailed state
CONDENSED_STATE = np.array([0, 2, 2, 4, 4, 1, 3, 5, 6, 7, 7, 6, 7, 7, 8])

#---------------------------------------------
# Methods
#---------------------------------------------
#---------------------------------------------
# Methods
#---------------------------------------------
def parse_command_line_args():
    '''Parse and validate command-line arguments.'''
    prog = os.path.basename(sys.argv[0])
    usage = 'Usage: %s [flags] <id_coef_file>\n\n' \
    'Calculate the ideal full and partial call rates of each sample.\n\n' \
    'Type ''%s -h'' to display full help.' % (prog, prog)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', action='store_true', dest='debug', default=False,
                      help='Print debugging information')
    (options, args) = parser.parse_args(sys.argv[1:])
    
    # Argument validation
    if len(args) != 1:
        print usage
        sys.exit(util.EXIT_BAD_INPUT_ARGS)
        
    # Read file name
    options.id_coef_file = args[0]
    return options

####################################################################################
class ImputationSimulation(object):
    '''Simulates gene dropping and measures call rates based on the selected identity states
    between a target sample ID i and each of the WGS training samples j in T.'''
    
    #---------------------------------------------
    # C-tors
    #---------------------------------------------
    def __init__(self, i, T, id_coef_file, e=0.01, debug=False, alpha=1.0):
        '''e=Desired relative error in call rate estimates.'''
        params = im.phase.PhaseParam(id_coef_file=id_coef_file)
        Delta = np.array([params.id_coefs(i, j)[1] for j in T])
        # p=row-stochastic matrix. Row j is the probability density of the detailed identity state
        # between i and j
        self.p = np.tile(SUB_PROB, (len(T), 1)) * Delta[:, CONDENSED_STATE]
        # Delta's might only approximately sum to 1, scale p to be a pdf
        self.p = statutil.scale_row_sums(self.p)
        self.debug = debug
        # Estimated # simulations for desired accuracy by the central-limit theorem
        self.num_simulations = int(np.ceil(1. / e ** 2))
        self.reset()
    #---------------------------------------------
    # Methods
    #---------------------------------------------
    def reset(self):
        '''Reset the simulation.'''
        # Full genotype call counter
        self.count_full = 0
        # Partial genotype call counter
        self.count_part = 0
        # Total # simulations
        self.count = 0

    def run(self, n=None):
        '''Run n simulations. If n=None, uses default #.'''
        n = n if n else self.num_simulations
        m = max(n / 10, 1000)
        for i in xrange(n):
            self._simulate_imputation()
            if self.debug and i % m == 0:
                print '\t%d/%d simulations completed ...' % (i, n)

    #---------------------------------------------
    # Properties
    #---------------------------------------------
    @property
    def call_rate_full(self):
        '''Full genotype call rate (%) estimate.'''
        return (1.0 * self.count_full) / self.count

    @property
    def call_rate_part(self):
        '''Partial genotype call rate (%) estimate.'''
        return (1.0 * self.count_part) / self.count
    
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------        
    def _simulate_imputation(self):
        '''Simulate imputing i from phased T-haplotypes.'''
        x = IMPUTABLE[statutil.multinomial_elementwise(self.p)]
        a0 = np.any(x / 2)  # Bit 0 set for at least one element of x ==> allele 0 is imputable
        a1 = np.any(x % 2)  # Bit 1 set for at least one element of x ==> allele 1 is imputable
        if a0 & a1:  # Both alleles are imputable
            self.count_full += 1
        if a0 | a1:  # Either allele is imputable
            self.count_part += 1
        self.count += 1

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    options = parse_command_line_args()
    try:
        T = im.examples.wgs_sample_id()
        I = im.examples.imputed_sample_id()
        index_of = im.examples.sample_id_to_index()
        n = len(I)
        # result = np.zeros((n, 4))
        e = 0.01
        for k, i in enumerate(I, 1):
            if options.debug:
                print 'Simulating sample %d (%d/%d)' % (i, k, len(I)) 
            sim = ImputationSimulation(i, T, options.id_coef_file, e=e, debug=options.debug)
            sim.run()
            # if options.debug:
                # print index_of[i], i, sim.call_rate_full, sim.call_rate_part
            # result[k] = [index_of[i], i, sim.call_rate_full, sim.call_rate_part]
            sys.stdout.write('%-4d %-10d %f %f\n' % (index_of[i], i, sim.call_rate_full, sim.call_rate_part))
        # np.savetxt(sys.stdout, result)
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(util.EXIT_FAILURE)
