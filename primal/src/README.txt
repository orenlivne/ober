====================================================================
Pedigree-Based Phasing and Imputation Software

Carole Ober and Dan Nicolae's Lab
Department of Human Genetics
The University of Chicago

Contact: Oren Livne <oren.livne@gmail.com>
====================================================================
This software performs phasing and imputation. This is sometimes called "Imputation Scenario A": given SNP data for all individuals and whole-genome-sequenced data at variants that includes the SNPs, impute all individuals at all variants. This can be broken into two separate steps 1 and 2 outlined below.

Requirements
------------
The code is written in Python. Version 2.7.3 (or a higher 2.x) required.

- The program runs separately for each autosomal chromosome 1..22. Sex chromosomes are not yet supported.
- It has been tested on Unix and Cygwin. 

***NOTE*** the test data is not available with this distribution because it takes 60M of space. It is stored in the Ober
Lab internal svn server (oberlab-tk).

- The parallel versions are suitable for the University of Chicago's Beagle Cray Supercomputer (http://beagle.ci.uchicago.edu/).
- You will need to install all necessarily python libraries. First, run system/bin/provision-python to install them in a separate virtual environment. Then, add the following lines to your .bashrc:

export APPS="<your chosen location for third-party applications, e.g., plink, tabix>"
export OBER="<your chosen location for a directory containing the source code, data and output dirs>"
source $OBER/system/dots/bash_profile

Individual directories (e.g. data and output) can also be overridden by exporting the appropriate variables before the "source" call. Refer to $OBER/system/dots/bash_profile for details.

1. Phasing
----------
input -> phase -> IBD segments -> index segments -> index_segments directory

Input:
file.{tped,tfam} - PLINK TPED SNP data on all individuals
file.kinship - kinship coefficients of all individual pairs
file.id - condensed identity coefficients of all individual pairs

Output:
- index_segments directory - index (for each SNP) of all IBD segments between all pairs of individuals

Serial programs:
impute/impute/phasing/phase.py
impute/bin/ibd_segments.py
impute/impute/ibd/index/segment_index.py

Parallel program:
impute/bin/pipeline-phase

2. Imputation
-------------
Input:
- index_segments directory (output of step 1)
- CGI data in tabixed tab-separated-values format on a subset of the individuals 

Output:
- Imputed data in tabixed tab-separated-values format on all individuals for all CGI variants

Serial programs:
impute/bin/run_chr_impute_cgi.py

Parallel program:
impute/bin/pipeline-phase
