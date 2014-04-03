'''
============================================================
Read and write Problem objects sets from/to a file set.
Support PLINK and NPZ formats, the former is a standard
in the genetics community, while the latter is much faster
(at least than loading/saving using nmp's loadxt
functions). 

Created on August 5, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import pickle, os, numpy as np, networkx as nx, util, db_gene
from impute.data import io_pedigree, io_genotype
from impute.data.problem import Problem
from impute.data.factory import GenotypeFactory
from impute.data.Pedigree import Pedigree
from db_gene.snp.ld_graph import Frames
from impute.data.constants import MISSING

#---------------------------------------------
# Methods
#---------------------------------------------
def read_plink(**kwargs):
    '''Load a problem from the following PLINK files:
    
        Default          Override Option    Data                                Format
        ======================================================================================
        prefix.pdg.tfam  pedigree           Pedigree adjacency                  PLINK TFAM
                                            (genotyped+nongenotyped samples)
        prefix.tfam      pedigree_genotyped Genotyped sample pedigree
                                            (sub-graph of the pedigree)         PLINK TFAM
                                            corresponding to prefix.tped
        prefix.tped      genotype           Genotype data                       PLINK TPED
        prefix.hap.tped  haplotype*         Haplotype data                      PLINK TPED 
        prefix.err       error**            Genotype errors flagged             Integer array (snps x samples) 
        prefix.info      info               Problem info                        pickle (binary)
        prefix.frm       frames             LD-independent SNP frames           text file
        prefix.lam       lam***             Haplotype est. recombination rate   text file
        
        * - hap data not loaded if this option is None.
        ** - errors set to 0 if this file is not found or this option is set to None.
        *** - data not loaded if if this file is not found.
    '''
    
    # Read input options
    verbose = kwargs.get('verbose', False)
    prefix = kwargs.get('prefix', '')
    overrideable_option = lambda name, default: kwargs.get(name, default if prefix else None)
    pedigree = overrideable_option('pedigree', prefix + '.pdg.tfam')
    pedigree_genotyped = overrideable_option('pedigree_genotyped', prefix + '.tfam')
    genotype = overrideable_option('genotype', prefix + '.tped')
    haplotype = overrideable_option('haplotype', prefix + '.hap.tped')
    error_file = overrideable_option('error', prefix + '.err')
    info = overrideable_option('info', prefix + '.info')
    if not np.all([[pedigree, pedigree_genotyped, genotype, error_file] is not None]):
        raise ValueError('Must specify a prefix or pedigree, pedigree_genotyped, genotype, error files')
    frames_file = overrideable_option('frames', prefix + '.frm')
    lam_file = overrideable_option('lam', prefix + '.lam')
    
    # Load data
    print_location = lambda x : x if x else '-'
    
    if verbose: print 'Reading pedigree from %s, %s ...' % (print_location(pedigree), print_location(pedigree_genotyped),)
    p = io_pedigree.read(pedigree, genotyped_id_file=pedigree_genotyped)
    
    if verbose:  print 'Reading genotype data from %s ...' % (print_location(genotype),)
    g = io_genotype.read('plink', 'genotype', tped=genotype, load_ids=False)
    
    if verbose: print 'Reading haplotype data from %s ...' % (print_location(haplotype),)
    h = io_genotype.read('plink', 'haplotype', tped=haplotype, load_ids=False) if haplotype else None
    
    if verbose: print 'Reading error data from %s ...' % (print_location(error_file),)
    error = np.loadtxt(error_file) if error_file and os.path.isfile(error_file) else None
    
    if verbose: print 'Reading frame data from %s ...' % (print_location(frames_file),)
    frames = db_gene.snp.ld_graph.read_frames(frames_file) if frames_file else None
    
    lam = np.loadtxt(lam_file) if lam_file and os.path.isfile(lam_file) else None
    
    # info = ProblemInfo(p, g) if info is None else info    
    problem = Problem(p, g, haplotype=h, error=error, frames=frames, lam=lam)
    if haplotype and info:
        if verbose: 
            print 'Reading problem info from %s ...' % (info,)
        with open(info, 'rb') as fout:
            problem.info = pickle.load(fout)
    return problem

def write_plink(problem, prefix, save_pedigree=True, save_node_type=True,
                save_genotype=True, save_haplotype=True, save_error=True, save_frames=True,
                verbose=False):
    '''Write problem to output file set in PLINK format:
    
        Default          Override Option    Data                                Format
        ======================================================================================
        prefix.pdg.tfam  pedigree           Pedigree adjacency                  PLINK TFAM
                                            (genotyped+nongenotyped samples)
        prefix.tfam      pedigree_genotyped Genotyped sample pedigree
                                            (sub-graph of the pedigree)         PLINK TFAM
        prefix.tped      genotype           Genotype data                       PLINK TPED
        prefix.hap.tped  haplotype          Haplotype data                      PLINK TPED 
        prefix.err       error              Genotype errors flagged             Integer array (snps x samples) 
        prefix.info      info               Problem info                        pickle (binary)
    '''
    if save_pedigree:
        if verbose: print 'Writing pedigree to %s.{tfam,.pdg.tfam} ...' % (prefix,)
        io_pedigree.write(problem.pedigree, prefix, save_node_type=save_node_type)

    if save_genotype:
        if verbose: print 'Writing genotype data to %s ...' % (prefix + '.tped',)
        with open(prefix + '.tped', 'wb') as fout: io_genotype.write('plink', problem.genotype, fout)

    if save_haplotype:
        if verbose: print 'Writing haplotype data to %s ...' % (prefix + '.hap.tped',)
        with open(prefix + '.hap.tped', 'wb') as fout: io_genotype.write('plink', problem.haplotype, fout)

    if save_error:
        if verbose: print 'Writing error data to %s ...' % (prefix + '.err',)
        with open(prefix + '.err', 'wb') as fout: np.savetxt(fout, problem.error, fmt='%d')

    if save_frames:
        if verbose: print 'Writing frames data to %s ...' % (prefix + '.frm',)
        with open(prefix + '.frm', 'wb') as fout: db_gene.snp.ld_graph.write_frames(problem.frames, fout)

    if verbose: print 'Writing problem info to %s ...' % (prefix + '.info',)
    with open(prefix + '.info', 'wb') as fout: pickle.dump(problem.info, fout)

####################################################################################
def read_npz(in_file):
    '''Read problem from NPZ file. in_file may be a file name or an open 
    file descriptor.'''

    files = np.load(in_file)
    graph = nx.DiGraph()
    graph.add_nodes_from(files['pedigree_nodes'])
    graph.add_edges_from(files['pedigree_graph'][0])
    p = Pedigree(graph,
                 sample_id=files['pedigree_sample_id'],
                 sex=files['pedigree_sex'],
                 phenotype=files['pedigree_phenotype'],
                 node_type=files['pedigree_node_type'],
                 sample_index=files['pedigree_sample_index'],
                 num_genotyped=files['pedigree_num_genotyped'][0])
    g = GenotypeFactory.new_instance('genotype', files['genotype_data'], files['genotype_snp'])
    h = GenotypeFactory.new_instance('haplotype', files['haplotype_data'], files['haplotype_snp'], qc=MISSING)
    error = files['error']
    h.qc = files['haplotype_qc']
    info = files['info'][0]
    frames = Frames((k, w) for k, v in files['frames'][0].iteritems() for w in v[0]) if files['frames'][0] else None
    lam = files['lam']
    
    # Optional fields
    if 'genotype_map' in files.files: g.map = files['genotype_map']
    if 'haplotype_poo_phase' in files.files: h.poo_phase = files['haplotype_poo_phase']
    if 'haplotype_hap_type' in files.files: h.hap_type = files['haplotype_hap_type']
    
    return Problem(p, g, haplotype=h, error=error, info=info, frames=frames, lam=lam)

def slim(problem):
    '''Remove large arrays from problem that are not necessary for downstream analysis: IBD
    segments (recalculated anyway with ibd_segments), QC, error.'''
    problem.haplotype.qc = []
    problem.error = []
    problem.info.ibd = []

def write_npz(problem, out_file):
    '''Write problem to NPZ file. out_file may be a file name or an open 
    file descriptor.'''
    p, g, h = problem.pedigree, problem.genotype, problem.haplotype
    if isinstance(out_file, str): util.mkdir_if_not_exists(os.path.dirname(out_file))
    # Wrap every non-np-array quantity by a np-array
    np.savez(out_file,
             
             pedigree_nodes=p.graph.nodes(),
             pedigree_graph=np.array([nx.to_edgelist(p.graph)]),
             pedigree_sample_id=p.sample_id,
             pedigree_sex=p.sex,
             pedigree_phenotype=p.phenotype,
             pedigree_node_type=p.node_type,
             pedigree_sample_index=p.sample_index,
             pedigree_num_genotyped=np.array([p.num_genotyped]),
             
             genotype_data=g.data,
             genotype_snp=g.snp,
             genotype_map=g.map,
             
             haplotype_data=h.data,
             haplotype_snp=h.snp,
             haplotype_qc=h.qc,
             haplotype_hap_type=h.hap_type,
             haplotype_poo_phase=h.poo_phase,
             
             error=problem.error,
             frames=np.array([problem.frames]),  # problem.frames.to_array(),
             info=np.array([problem.info]),
             lam=problem.lam)

def plink_to_npz(prefix, npz_file, pedigree_file=None, pedigree_genotyped=None, frames=None,
                 verbose=False, **kwargs):
    '''Convert Problem plink -> npz.'''
    if verbose: print 'Reading plink input...'
    haplotype = kwargs.get('haplotype', prefix + '.hap.tped')
    problem = read_plink(pedigree=pedigree_file if pedigree_file else prefix + '.pdg.tfam',
                         pedigree_genotyped=pedigree_genotyped if pedigree_genotyped else prefix + '.tfam',
                         genotype=prefix + '.tped',
                         haplotype=haplotype,
                         error=prefix + '.err',
                         info=prefix + '.info',
                         frames=frames if frames else prefix + '.frm',
                         lam=prefix + '.lam',
                         verbose=verbose)
    if verbose: print 'Writing npz output to %s...' % (npz_file,)
    write_npz(problem, npz_file)
    return problem
    
def npz_to_plink(npz_file, prefix, verbose=False):
    '''Convert Problem npz -> plink.'''
    if verbose: print 'Reading npz input...'
    problem = read_npz(npz_file)
    if verbose: print 'Writing plink output...'
    write_plink(problem, prefix, verbose=verbose)
    return problem

def write_family_txt(p, out_file, f=None, delimiter=' '):
    '''Save family f in problem p to file name out_file in text format. Row format:
    father_genotype 
    mother_genotype
    father_haplotype
    mother_haplotype
    child1_haplotype
    ...
    childN_haplotype
    If f is None, uses the first genotyped family.'''
    if f is None: f = p.families(genotyped=True)[0]
    np.savetxt(out_file, np.concatenate((np.reshape(p.genotype.data [:, f.parents, :],
                                                    (p.num_snps, 4)),
                                         np.reshape(p.haplotype.data[:, f.parents, :],
                                                    (p.num_snps, 4)),
                                         np.reshape(p.haplotype.data[:, f.children_list, :],
                                                    (p.num_snps, 2 * len(f.children_list)))),
                                        axis=1), fmt='%d', delimiter=delimiter)

def read_info_npz(in_file):
    '''Read a ProblemInfo from NPZ file. in_file may be a file name or an open file descriptor.'''
    return np.load(in_file)['info'][0]
    
def write_info_npz(info, out_file):
    '''Write a ProblemInfo object into an NPZ file.'''
    np.savez(out_file, info=np.array([info]))

#---------------------------------------------
# Private Methods
#---------------------------------------------
