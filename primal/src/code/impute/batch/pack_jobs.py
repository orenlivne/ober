#!/usr/bin/env python
'''
============================================================
Creates CRI/Beagle PBS job submission scripts to span a
parameter space.

Created on February 4, 2013
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import sys, os, itertools as it, traceback, ConfigParser, optparse, re, stat, subprocess as sub

#---------------------------------------------
# Constants
#---------------------------------------------
# This program's name
PROGRAM = os.path.basename(sys.argv[0])
# No. cores per Beagle node
CORES_PER_NODE = 24

'''Abnormal main program termination codes'''
EXIT_FAILURE = 1
EXIT_BAD_INPUT_ARGS = 2
EXIT_FILE_NOT_FOUND = 3

#---------------------------------------------
# Methods
#---------------------------------------------
def mkdir_if_not_exists(path):
    '''Create a directory if it does not exist yet.'''
    # print 'Making directory', path
    if path and not os.path.exists(path):
        os.makedirs(path)
        
def run_command(cmd, die_on_exception=True, verbose=False):
    '''Run command in a sub-shell. Return the command's return code.'''
    if verbose:
        print 'Running', cmd
    # subprocess.call(cmd)
    p = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
    output, errors = p.communicate()
    if p.returncode != 0:
        if verbose:
            print 'Output', output
            print 'Errors', errors
        if die_on_exception:
            raise ValueError('Command failed:', cmd)
    return p.returncode, output, errors

class mdict(dict):
    '''Multivalued dictionary from 
    http://code.activestate.com/recipes/440502-a-dictionary-with-multiple-values-for-each-key/'''
    def __setitem__(self, key, value):
        '''Add the given value to the list of values for this key'''
        self.setdefault(key, []).append(value)
        
    @staticmethod
    def from_items(items):
        '''Initialize from a key-value tuple list.'''
        a = mdict()
        for k, v in items:
            a[k] = v
        return a

def cumsum(it):
    '''Cumulative sum. Available in numpy, but we cannot rely on it being available.'''
    total = 0
    for x in it:
        total += x
        yield total
        
def prod(it):
    '''A product of the iterates. Available in numpy, but we cannot rely on it being available.'''
    return reduce(lambda x, y: x * y, it)

def range_sizes(n, p):
    '''Break the integers 0..n-1 into p more-or-less equal ranges. This is an iterator over ranges.
    Usage:
    >>> range_sizes(10,3) = [4,3,3]'''
    q, r = n / p, n % p
    return it.chain(it.islice(it.repeat(q + 1), r), it.islice(it.repeat(q), p - r))

def parse_param(name, value):
    '''Parses a configuration value string into the proper Param sub-class.
    TODO: replace by a more sophisticated and flexible parser.'''
    value_regex = re.compile('\s*(\d+)\s*:\s*(\d+)\s*')
    m = value_regex.match(value)
    if m:
        # Range, e.g., 123:456. Inclusive on both ends.
        return DiscreteParam(name, range(int(m.group(1)), int(m.group(2))))
    else:
        # Discrete range
        values = re.compile('\s*,\s*').split(value.strip(' '))
        return DiscreteParam(name, values)

def pbs_script_writer(target, config, out_file):
    '''A factory method of PbsScriptWriter objects.'''
    if target == 'beagle':
        return BeagleScriptWriter(config, out_file)
    elif target == 'cri':
        return JobArrayScriptWriter(config, out_file)
    else:
        raise ValueError('Unsupported target system ''%s''' % (target,))

####################################################################################
class ParamIterator(object):
    '''Iterates over a cross-product of parameter ranges. Outputs parameter combinations
    in blocks so that there are total of p blocks (e.g., p=#cores in a parallel system).'''

    def __init__(self, params, p, predicate=None):
        '''Initialize from a list of parameter ranges.'''
        self.params = params
        self.p = p
        self.params = params
        # Calculate total # of combinations by looping over all parameters
        self.size = len(list(self._generate_param_iter(params, predicate)))
        # Store iterator for subsequent iteration
        self.param_iter = self._generate_param_iter(params, predicate)
        self.block_sizes = range_sizes(self.size, self.p)

    def __iter__(self):
        return self
    
    @staticmethod
    def _generate_param_iter(params, predicate):
        param_iter = it.product(*params)
        if predicate:
            if len(params) == 2:
                param_iter = it.ifilter(ParamIterator._predicate_filter(predicate), param_iter)
            else:
                raise ValueError('Predicate supports exactly two parameters')
        return param_iter
    
    @staticmethod
    def _predicate_filter(predicate):
        '''A filter corresponding to a predicate. Must have exactly two parameters for this to work.'''
        # TODO: replace with a full-fledged parser
        if predicate == '<':
            return lambda x: x[0] < x[1]
        elif predicate == '<=':
            return lambda x: x[0] <= x[1]
        elif predicate == '>':
            return lambda x: x[0] > x[1]
        elif predicate == '>=':
            return lambda x: x[0] >= x[1]
        else:
            raise ValueError('Unsupported predicate ''%s''' % (predicate,))

    def next(self):  # @ReservedAssignment
        '''An iterator of parameter combinations.'''
        for block_size in self.block_sizes:
            return [self.param_iter.next() for _ in xrange(block_size)]
        raise StopIteration

####################################################################################
def group_iterates(it, c):
    '''Groups ParamIterator iterator (it) results into blocks of size c each
    (e.g., c=#cores per node on Beagle). Iterates are also enumerated within the block.'''
    buf = []
    for k, x in enumerate(it):
        buf.append((k, x))
        if k % c == (c - 1):
            yield buf
            buf = []
    if buf:
        yield buf

####################################################################################
class Param(object):
    '''Represents an input parameter and its range.'''

    def __init__(self, name):
        '''Initialize from a comma-delimited string.'''
        self.name = name
    
    def __repr__(self):
        return '%s[%s]' % (self.__class__, self.name,)

    @property
    def size(self):
        '''Total number of parameter values.'''
        raise ValueError('Must be implemented by sub-classes.')
    
    def __iter__(self):
        '''An iterator over the parameter''s range of values.'''
        raise ValueError('Must be implemented by sub-classes.')

####################################################################################
class DiscreteParam(Param):
    '''Represents an input parameter ranging over a list of values.'''

    def __init__(self, name, values):
        '''Initialize from a list of values.'''
        Param.__init__(self, name)
        self.name = name
        self.values = values
    
    @property
    def size(self):
        '''Total number of parameter values.'''
        return len(self.values)
    
    def __iter__(self):
        '''An iterator over the parameter''s range of values.'''
        return iter(self.values)
     
####################################################################################
class Configuration(object):
    '''Configuration file parser.'''
    
    def __init__(self, config, options=None):
        # Add reserved variables
        self._config = config
        self.options = options
        
        # Init parameter list, cached properties
        self.instances_per_node = int(self['pbs']['instances_per_node'])
        self.processes = int(self['pbs']['processes'])
        self.nodes = int(self['pbs']['nodes'])
        self.params = [parse_param(k, v) for k, v in self['param'].iteritems() if not k.startswith('_')]
        if self.params:
            print 'Recognized %d parameter(s): %s' % (len(self.params), ', '.join(param.name for param in self.params))
        # Finish initialization        
        self.validate()

    @staticmethod
    def from_file(file_name, options=None, params=None):
        '''A factory method. Parse a file into a Configuration object.
        The options dictionary is added in the ''options'' section.'''
        config = ConfigParser.ConfigParser()
        config.read(config_file)
        
        # Add options
        config.add_section('options')
        for k, v in ((k, v) for k, v in options.__dict__.iteritems() if v):
            config.set('options', k, v)

        # Add parameters
        if params:
            for k, v in params.iteritems():
                config.set('DEFAULT', '_' + k, v)

        return Configuration(config, options=options)

    def __iter__(self):
        if self['param'].has_key('_executable'):
            # External parameter splitting command specified, create a 1-D parameter space
            # (one parameter per core)
            return ParamIterator([parse_param('index', '0:%d' % (self.cores))], self.cores)
        else:
            try:
                predicate = self._config.get('param', '_predicate')
            except ConfigParser.NoOptionError:
                # No parameters
                predicate = None
            return ParamIterator(self.params, self.cores, predicate=predicate)
    
    def __getitem__(self, section):
        '''Convert a config section to a dictionary.
        @see http://wiki.python.org/moin/ConfigParserExamples'''
        dict1 = {}
        options = self._config.options(section)
        for option in options:
            dict1[option] = self._config.get(section, option)
        return dict1
    
    @property
    def defaults(self):
        return self._config.defaults()
    
    @property
    def cores(self):
        return self.nodes * self.instances_per_node
    
    @property
    def modules(self):
        return re.compile('\s*,\s*').split(self['env']['modules'].strip(' ')) 
    
    @property
    def job_name(self):
        return self['pbs']['job_name']
    
    @property
    def input_files(self):
        '''A dictionary of input file-to-be-copied - to - its-alias-under-the-mem-directory.
        Read from a comma-delimited string such as ''a:b, c, d:e''.
        If an item not in the format a:b, its alias is set to its base name.'''
        return dict((parts[0], parts[-1] if len(parts) == 2 else os.path.basename(parts[0])) 
                    for parts in (x.split(':') for x in re.compile('\s*,\s*').split(self['exec']['transfer_input_files'].strip(' ')) if x))
    
    @property
    def output_files(self):
        return [x for x in re.compile('\s*,\s*').split(self['exec']['transfer_output_files'].strip(' ')) if x]
    
    @property
    def jobs(self):
        '''Return the number of jobs = # parameter configurations.'''
        return self.__iter__().size

    def validate(self):
        '''Input validation.'''
        i = self.instances_per_node
        if i < 1 or i > CORES_PER_NODE:
            print usage
            print('\nMust specify # jobs per node in in 1..%d. Found %s' % (CORES_PER_NODE, i))
            sys.exit(EXIT_BAD_INPUT_ARGS)
        jobs = self.jobs
        if self.cores > jobs:
            print '#jobs %d, #cores %d, #cores/node %d, Adjusting # nodes' % (jobs, self.cores, i)
            self.nodes = (jobs - 1) / i + 1
            if self.nodes == 1 and self.cores > jobs:
                print 'Adjusting # cores as well'
                self.instances_per_node = self.jobs
        print 'Using %d nodes, %d cores/node, %d cores' % (self.nodes, self.instances_per_node, self.cores)
        print '# #Jobs %d, jobs/core %d' % (jobs, (jobs - 1) / self.cores + 1)

####################################################################################   
class FileWriter(object):
    '''An abstract class that writes an output file.'''

    def __init__(self, out_file):
        self.out_file = out_file
        
    def write_file(self):
        '''Write the file.'''
        raise ValueError('Sub-classes must implement this method')
    
    def write(self):
        '''A a template method - main write call.'''
        with open(self.out_file, 'wb') as self.out:
            self.write_file()

    def write_header(self):
        '''Write a file header.'''
        # Header
        self.writeln('#' * 75)
        self.writeln('# %s: file generated by %s' % (os.path.basename(self.out_file), PROGRAM))
        self.writeln('#' * 75)
        self.writeln('')

    def writeln(self, s):
        '''Write a string to file and terminate with a new line.'''                    
        self.out.write(s + '\n')

####################################################################################   
class BeagleNodeFileWriter(FileWriter):
    '''Writes a job file to be run on a single node. This is a bash script that packs copies of the
    executable specified in the config file'''
    
    def __init__(self, out_file, config, node, group, write_in_files):
        FileWriter.__init__(self, out_file)
        # Create a defensive copy of the configuration
        self.config = config
        self.node = node
        self.group = group
        self.out_dir = os.path.dirname(self.out_file)
        self.write_in_files = write_in_files
        mkdir_if_not_exists(self.out_dir)

    def write_file(self):
        '''Write the node job file. This is a bash script that packs copies of the
        executable specified in the config file.'''
        
        # Inject parameter combination representation into configuration for string interpolation
        node_name = '%04d' % (self.node,)
        self.config._config.set('DEFAULT', '_node', node_name)
        self.config._config.set('DEFAULT', '_param', 'None')

        # Set useful bash script alias variables, inject into string interpolation
        self.writeln('#!/bin/bash')
        self.write_header()
        try: _mem = self.config.defaults['_mem']
        except KeyError: _mem = self.out_dir
        self.writeln('ramdisk="%s"' % (_mem,))
        ramdisk = '${ramdisk}'
        self.config._config.set('DEFAULT', '_mem', ramdisk)
        
        # Input file copying section
        self.writeln('# Clear RAMDISK memory')
        # self.writeln('transfer_input_files %s' % (self.config['exec']['transfer_input_files'],))
        # self.writeln('%s' % (repr(self.config.input_files),))
        self.writeln('if [ ! -d %s ]; then' % (ramdisk,))
        self.writeln('  mkdir -p %s' % (ramdisk,))
        self.writeln('fi')
        if self.config.input_files:
            self.writeln('find %s -mindepth 1 -maxdepth 1 -user $USER -exec rm -rf {} \;' % (ramdisk,))
        self.writeln('# Copy input files to RAMDISK memory')
        # Copy input files to their alias locations under the ramdisk directory
        for k, v in self.config.input_files.iteritems():
            parent_dir = os.path.dirname(v)
            if parent_dir:
                alias = ramdisk + '/' + v
                self.writeln('mkdir -p %s' % (os.path.dirname(alias),))
            self.writeln('cp -r %s %s/%s' % (k, ramdisk, v))
            # Make sure create directory is writable so that we can remove it later
#            if parent_dir:
#                self.writeln('chmod -R +w %s' % (alias,))
        self.writeln('')
        
        self.writeln('# Spawn jobs')
        for k, params in self.group:
            # Inject parameter combination representation into configuration for string interpolation
            params_name = '%04d' % (k,)
            self.config._config.set('DEFAULT', '_param', params_name)
            in_file_name = '%s/%s-%s.in' % (self.out_dir, self.config.job_name, params_name)
            if self.write_in_files:
                with open(in_file_name, 'wb') as in_file:
                    for param in params:
                        in_file.write(' '.join(map(lambda x: '%s' % (str(x),), param)) + '\n')
            if self.config.options.manual_input:
                self.writeln('%s &' % (self.config['exec']['executable'],))
            else:
                self.writeln('%s 2> %s/%s-%s.log < %s &' % \
                             (self.config['exec']['executable'],
                              self.out_dir, self.config.job_name, params_name, in_file_name))
#            self.writeln('%s 2> /dev/null < %s &' % \
#                         (self.config['exec']['executable'], in_file_name))
        # Add a CPU and Memory monitoring thread
#         self.writeln('( top -b -d%d -u oren | tee >(grep Mem) >(grep Cpu) >& /dev/null ) >& %s/%s.monitor &' % \
#                      (self.config.options.frequency, self.out_dir, self.config.job_name))
        if self.config.options.frequency > 0:
            self.writeln('monitor-usage %d false >& %s/%s.monitor &' % (self.config.options.frequency, self.out_dir, self.config.job_name))
        self.writeln('wait')
        self.writeln('')
        
        self.writeln('# Clean RAMDISK memory')
        if self.config.output_files: 
            self.writeln('mv %s %s' % (' '.join(map(lambda x : '%s/%s' % (ramdisk, x), self.config.output_files)),
                                       self.out_dir))
        # Delete all top-level directories created under ramdisk
        if self.config.input_files:
            self.writeln('find %s -mindepth 1 -maxdepth 1 -user $USER -exec rm -rf {} \;' % (ramdisk,))
        
        # Restore config variables set for node string interpolation
        self.config._config.set('DEFAULT', '_mem', _mem)

####################################################################################   
class PbsScriptWriter(FileWriter):
    '''Writes the main Beagle PBS script for all IBD segment jobs to an output stream.'''

    def __init__(self, config, out_file):
        '''Prepare a writer that will submit jobs to ''nodes''-Beagle nodes.'''
        FileWriter.__init__(self, out_file)
        self.config = config
        self.c = self.config.instances_per_node
        self.jobs = self.config.jobs
     
    def _custom_headers(self):
        '''A dictionary of custom PBS headers to write.'''
        pass

    def _write_environment(self):
        '''Initialize the execution environment.'''
        pass
    
    def write_body(self):
        '''Write the body of the PBS submission script.'''
        raise ValueError('Must be implemented by sub-classes')

    def write_file(self):
        '''Main call that writes the PBS script to the file ''out_file''.'''
        cores = self.config.cores
        self.writeln('# Using %d nodes, %d cores/node, %d cores' % (self.config.nodes, self.c, cores))
        self.writeln('#Jobs %d, jobs/core %d' % (self.jobs, (self.jobs - 1) / cores + 1))

        # PBS headers
        self.writeln('#!/bin/bash')
        for k, vals in mdict.from_items(self._common_headers() + self._custom_headers()).iteritems():
            for v in vals:
                self.writeln('#PBS %s %s' % (k, v))
        self.writeln('')
        
        self._write_environment()
        self.writeln('cd $PBS_O_WORKDIR')
        self.writeln('')
        
        self.write_body()

    def _common_headers(self):
        headers = [
            ('-N', self.config.job_name),
            ('-j', 'oe'),
            ('-l', 'walltime=%s' % (self.config['pbs']['walltime'],)),
#                ('-l', 'nodes=1'),
#                ('-l', 'ppn=%d' % (config.instances_per_node,))
            ]
        queue = self.config['pbs']['queue']
        if queue:
            headers += [('-q', self.config['pbs']['queue'])]
        return headers

####################################################################################   
class BeagleScriptWriter(PbsScriptWriter):
    '''Writes the main Beagle PBS script for all IBD segment jobs to an output stream.'''

    def __init__(self, config, out_file):
        '''Prepare a writer that will submit jobs to ''nodes''-Beagle nodes.'''
        PbsScriptWriter.__init__(self, config, out_file)

    def _custom_headers(self):
        '''A dictionary of custom PBS headers to write.'''
        return [
                ('-l', 'mppwidth=%d' % (CORES_PER_NODE * self.config.nodes,)),
                ('-A', self.config['pbs']['project'])
                ]

    def _write_environment(self):
        '''Initialize the module script.'''
        self.writeln('echo /opt/modules/default\n' \
            '. /opt/modules/default/init/bash\n' \
            'module swap PrgEnv-pgi PrgEnv-gnu')
        for module in self.config.modules:
            self.writeln('module load %s' % (module,))
        self.writeln('module list 2>&1')

    def write_body(self):
        '''Write the body of the PBS submission script.'''
        split_executable = self.config['param'].has_key('_executable')
        if split_executable:
            # Run split command that is supposed to generate in files. Do not write in files below.
            split_cmd = self.config['param']['_executable']
            if split_cmd == 'None':
                print 'Skipping split command'
            else:
                print 'Running split command: %s' % (split_cmd,)            
                run_command(split_cmd)
            
        for node, group in enumerate(group_iterates(self.config, self.c)):
            # print 'node', node, 'cores', group[0][0], 'to', group[-1][0]
            node_file = self._node_file_name(node, group)
            # Assuming executable is on /lustre; bypassing executable transfer to compute nodes
            self.writeln('aprun -b -n 1 -N 1 -d %d %s &' % (self.config.processes, node_file))
            BeagleNodeFileWriter(node_file, config, node, group, not split_executable).write()
            os.chmod(node_file, stat.S_IRWXU)
        self.writeln('wait')

    def _node_file_name(self, node, group):
        '''Return the node script file name.'''
        out_dir = os.path.dirname(self.out_file)
        return '%s/node-%04d/%s-%04d_to_%04d.sh' % (out_dir, node, self.config.job_name, group[0][0], group[-1][0])
        
####################################################################################   
class JobArrayScriptWriter(PbsScriptWriter):
    '''Writes the main CRI Torque job array script for all IBD segment jobs to an output stream.'''

    def __init__(self, config, out_file):
        '''Prepare a writer that will submit jobs to ''nodes''-Beagle nodes.'''
        config.instances_per_node = 1
        PbsScriptWriter.__init__(self, config, out_file)
        
    def _custom_headers(self):
        '''A dictionary of custom PBS headers to write.'''
        return [('-t', '0-%d' % (self.config.nodes - 1,))]
        
    def _write_environment(self):
        '''Initialize the module script.'''
        self.writeln('n=${PBS_ARRAYID} # A convenient alias to the node number')

    def write_body(self):
        '''Write the body of the PBS submission script.'''
        # Body is parameterized by the PBS_ARRAYID parameter, i.e., the node number. So only
        # one section is required for all nodes, not one per node as in the Beagle implementation.
        mem = self.config.defaults['_mem']
        node_name = '$n'
        ramdisk = '%s/%s-node-%s' % (mem, self.config.job_name, node_name)  # Node ramdisk directory
        self.config._config.set('DEFAULT', '_mem', ramdisk)
        self.config._config.set('DEFAULT', '_node', node_name)
        self.config._config.set('DEFAULT', '_param', 'None')
        # Create a temp dir under RAMDISK that is specific to our job. TODO: name this uniquely to prevent conflicts?

        self.writeln('# Copy input files to RAMDISK memory')
        self.writeln('mkdir -p %s' % (ramdisk,))
        self.writeln('cp -r %s %s' % (' '.join(self.config.input_files), ramdisk))
        self.writeln('')
        
        node_dir = self._node_dir(node_name)
        in_file_name = self._in_file_name(node_name)
        self.writeln('# Spawn job')
        self.writeln('%s 2> %s/%s-%s.log < %s' % \
                     (self.config['exec']['executable'],
                      node_dir, self.config.job_name, node_name, in_file_name))
        self.writeln('')

        # Creaet input files        
        for node, group in enumerate(group_iterates(self.config, self.c)):
            # print 'node', node
            node_name = '%d' % (node,)  # Must match job array ID so that we can use it here
            self._write_node_inpute_file(node, group, self._in_file_name(node_name))
                
        self.writeln('# Clean RAMDISK memory')
        if self.config.output_files: 
            self.writeln('mv %s %s' % (' '.join(map(lambda x : '%s/%s' % (ramdisk, x), self.config.output_files)),
                                        out_dir))
        self.writeln('rm -rf %s' % (ramdisk,))
        self.writeln('')
        
        # Restore state
        self.config._config.set('DEFAULT', '_mem', mem)

    def _write_node_inpute_file(self, node, group, in_file_name):
        '''Write job input file.'''
        mkdir_if_not_exists(os.path.dirname(in_file_name))
        _, params = list(group)[0]
        with open(in_file_name, 'wb') as in_file:
            for param in params:
                in_file.write(' '.join(map(lambda x: '%s' % (str(x),), param)) + '\n')

    def _node_dir(self, node_name):
        out_dir = os.path.dirname(self.out_file)
        return '%s/node-%s' % (out_dir, node_name)

    def _in_file_name(self, node_name):
        return '%s/%s-node-%s.in' % (self._node_dir(node_name), self.config.job_name, node_name)

####################################################################################
if __name__ == '__main__':
    '''
    --------------------------------------------------
    Main program
    --------------------------------------------------
    '''
    # Parse command-line arguments
    usage = 'Usage: %s [flags] <config-file> <out-dir>\n\n' \
        'Create CRI cluster or Beagle PBS job submission scripts for spanning a\n' \
        'parameter space.\n\n' \
        'Config parameters (''-p'' flag) become global placeholders in the config file whose\n' \
        'keys are prefixed with ''_''. For instance, -p chrom=22 results in a placeholder\n' \
        '_chrom whose value is 22.\n\n' \
        'Type ''%s -h'' to display full help.' % (PROGRAM, PROGRAM)
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-v', '--debug'          , action='store_true'  , dest='debug', default=False,
                      help='Print debugging information')
    parser.add_option('-t', '--target', type=str, dest='target', default='cri',
                      help='Target system (beagle|cri)')
    parser.add_option('-p', '--params', type=str, dest='params', default='',
                      help='Config parameters in the format key1=value1,key2=value,...,keyN=valueN')
    parser.add_option('-f', '--monitor-frequency', type='int', dest='frequency', default=0,
                      help='CPU & memory thread monitoring frequency [sec]. If 0, is turned off')
    parser.add_option('-m', '--manual-input'          , action='store_true'  , dest='manual_input', default=False,
                      help='Specify input stream manually in sub file instead of generating it')
    options, args = parser.parse_args(sys.argv[1:])
    # print options
    # print args
    if len(args) != 2:
        print usage
        sys.exit(EXIT_BAD_INPUT_ARGS)
    config_file, out_dir = args[0], os.path.abspath(args[1])
    if not os.path.isfile(config_file):
        print 'Config file %s not found' % (config_file,)
        sys.exit(EXIT_FILE_NOT_FOUND)
    # Parse config file
    params = dict(x.split('=') for x in options.params.split(',')) if options.params else {}
    config = Configuration.from_file(config_file, options, params)
    
    # Main program that generates the PBS script
    try:
        mkdir_if_not_exists(out_dir)
        writer = pbs_script_writer(options.target, config, out_dir + '/' + config.job_name + '.pbs')
        writer.write()
    except:
        traceback.print_exc(file=sys.stdout)
        sys.exit(EXIT_FAILURE)
