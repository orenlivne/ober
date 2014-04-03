'''
============================================================
Phasing algorithm. Implements the chain-of-responsibility
pattern where filter=phasing stage (or phaser, for short)
operating on request=a Haplotype object.

This class contains the main processing chain. See other
phase_* modules for specific parts of the chain.
 
Created on July 26, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from impute.phasing.phase_trivial import trivial_phaser
from impute.phasing.phase_family import family_child_comparison_phaser, family_phaser
from impute.phasing.phase_core import new_phaser_chain
from impute.phasing.pre_processing import prepare_phaser
from impute.phasing.phase_distant import family_sib_comparison_phaser

'''Main phasing processing chain within a family.'''
def main_phaser(debug=False, print_times=False):
    return new_phaser_chain([
                             prepare_phaser(),
                             trivial_phaser(),
                             family_phaser(),
                             family_child_comparison_phaser(),
                             family_sib_comparison_phaser(),
                             ], debug=debug, print_times=print_times)

'''Only phase within nuclear families.''' 
def nuclear_family_phaser(debug=False, print_times=False):
    return new_phaser_chain([
                             prepare_phaser(),
                             trivial_phaser(),
                             family_phaser(),
                             family_child_comparison_phaser(),
                             family_sib_comparison_phaser(),
                             ], debug=debug, print_times=print_times)

def phaser_stages13(debug=False, print_times=False):
    return new_phaser_chain([
                             prepare_phaser(),
                             trivial_phaser(),
                             family_phaser(),
                             family_child_comparison_phaser(),
                             ], debug=debug, print_times=print_times)

def prepare_and_trivial_phaser(debug=False, print_times=False):
    return new_phaser_chain([
                             prepare_phaser(),
                             trivial_phaser()
                             ], debug=debug, print_times=print_times)
