#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
"""
WholeGenomePhylogeny.py is a pipeline for maximum likelihood phylogenetic inference from fasta formated genomes. It is written to run on the Berkeley Research Computing's (BRC) High Performance Computing (HPC) cluster, using the SLURM scheduler for process management and queueing within the shared computing environment.

Use option -h/--help for usage information.

Dependencies:
 - Python 2.7
 - Biopython
 - NumPy
 - MUMmer
 - MAFFT
 - PartitionFinder2
 - RAxML
 - ExaML
 - OpenMPI
 - ETE toolkit version 3

Using the Anaconda2 distribution of Python is highly recommended. PartitionFinder2, in particular, relies on many non-standard Python modules. Conflicting installs of Python can cause errors when importing modules, so within the BRC HPC use 'module purge' prior to running.

Requires the following modules on the BRC HPC: raxml/8.1.17 gcc/4.8.5 openmpi/1.10.2-gcc
"""
# Standard modules
import argparse, os, sys, time, pickle
import subprocess as sp

# Nonstandard modules

# My modules
import genome_datasets
from WGP1_orthology import orthology
from WGP2_multiple_alignment import multiple_alignment
from WGP3_partition import partition
from WGP4_phylogeny import phylogeny
from WGP5_bootstrap import bootstrap

__author__ = "Christopher Hann-Soden"
__copyright__ = "Copyright 2016, Christopher Hann-Soden"
__credits__ = ["Christopher Hann-Soden"]
__licence__ = "GPL"
__version__ = "0.8"
__maintainer__ = "Christopher Hann-Soden"
__email__ = "channsoden@berkeley.edu"
__status__ = "Development"

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generates a whole genome maximum likelihood phylogeny from a list of fasta genomes.')
    parser.add_argument('-r', '--reference', type=str, required=True, help='path to reference genome')
    parser.add_argument('-s', '--step', type=int, default=1, help='step of pipeline to start on')
    parser.add_argument('-e', '--seed', type=int, default=time.time(), help='set pseudo-random seed for reproducable results')
    parser.add_argument('-o', '--output', type=str, default='WGP_analysis', help='name stub for output files')
    parser.add_argument('-d', '--dataset', type=str, default='', help='name of a dictionary of fasta genomes in genome_datasets.py')
    parser.add_argument('-b', '--bootstrap', dest='bootstrap', action='store_true', help='run a bootstrap analysis after making tree')
    parser.add_argument('genomes', type=str, nargs='*', help='list of fasta formatted genomes')
    parser.set_defaults(bootstrap=False)
    args = parser.parse_args()

    if args.genomes and args.dataset:
        raise Exception('Either provide a list of fasta genomes or name a dataset with -d, not both.')
    elif not args.genomes and args.dataset:
        try:
            import genome_datasets
            args.genomes = eval('genome_datasets.'+args.dataset).values()
        except ImportError:
            raise ImportError('genome_datasets.py was not found in your PYTHONPATH')
        except AttributeError:
            raise AttributeError(args.dataset+' not found in genome_datasets.py')
    elif args.genomes:
        pass
    else:
        raise Exception('Either a dataset must be provided with -d or a list of fasta formatted genomes.')
        
    return args

def main(args):
    with open(args.output+'.args.pickle', 'wb') as fh:
        pickle.dump(args, fh)

    if args.step == 1:
        print 'Began orthology search'
        seg_files = orthology(args)
        if not seg_files:
            # Orthology search must have failed
            return None
        print 'Finished orthology search'
        os.remove('1_orthology/'+args.output+'.uni_shared_ref.pickle')
        args.step += 1

    if args.step == 2:
        print 'Began multiple alignment'
        alignment = multiple_alignment(args, seg_files)
        print 'Finished multiple alignment'
        args.step += 1

    if args.step == 3:
        print 'Began partitioning'
        partition_file, phylip = partition(args, alignment)
        print 'Finished partitioning'
        args.step += 1

    if args.step == 4:
        print 'Began ML tree optimization'
        tree = phylogeny(args, partition_file, phylip)
        print 'Finished ML tree optimization'
        args.step += 1

    if args.step == 5:
        if args.bootstrap:
            bs_tree = bootstrap(args, alignment, tree)
        else:
            bs_tree = tree
        args.step += 1
        
    print 'Completed whole genome phylogenetic analysis.'
    os.remove(args.output+'.args.pickle')
    return bs_tree

if __name__ == '__main__':
    args = parse_arguments()
    main(args)
