#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import os
from argparse import Namespace

# Nonstandard modules
from ete3 import Tree
import numpy as np

# My modules
from processing_tools import mapPool
from SLURM_tools import submit, job_wait
from fasta_tools import fasta_to_dict
from WGP3_partition import partition
from WGP4_phylogeny import phylogeny
from WGP6_cleanup import cleanup

def bootstrap(args, alignment, tree, reps = 100):
    bs_args = range_args(args, reps)
    bs_alignments = sample(alignment, reps)

    partition_jobs = [(partition, (args, bsa) ) for args, bsa in zip(bs_args, bs_alignments)]
    part_results = mapPool(100, partition_jobs)

    phylo_jobs = [(phylogeny, (args, pr[0], pr[1])) for args, pr in zip(bs_args, part_results)]
    ID = submit(phylo_jobs,
                pool=True,
                job_name='ExaML',
                modules=['raxml/8.1.17'])
    job_wait(ID)
    bs_trees = get_trees(bs_args)

    bs_tree = map_support(tree, bs_trees)

    [cleanup(args) for args in bs_args]

    return bs_tree

def range_args(args, reps):
    new_args = []
    for i in range(reps):
        new = Namespace(**vars(args))
        args.output = args.output+'.BS'+str(i)
        new_args.append(new)
    return new_args

def sample(alignment, reps):
    records = fasta_to_dict(alignment)
    names, seqs = zip(*records.items())
    aln = np.array([list(seq) for seq in seqs])
    M = aln.shape[1]

    name_lens = [len(name) for name in names]
    longest = max(name_lens)
    heads = np.array([['>']+list(name)+[' ']*(longest-l)+['\n'] for name, l in zip(names, name_lens)])
    tails = np.array([['\n'] for name in names])
    
    outfiles = []
    for i in range(reps):
        sample_idxs = np.random.choice(M, size=M, replace=True)
        new_aln = aln.T[sample_idxs].T
        filename = '{}.BS{}.fa'.format(alignment.rsplit('.', 1)[0], i)
        with open(filename, 'w') as fh:
            fh.write( np.concatenate((heads, new_aln, tails), axis=1).tostring() )
        outfiles.append(filename)

    return outfiles
    
def get_trees(bs_args):
    trees = ['4_phylogeny/'+'ExaML_result.'+args.output for args in bs_args]

def map_support(tree_file, bs_tree_files):
    tree = Tree(tree_file)
    bs_trees = [Tree(bst) for bst in bs_tree_files]
    
    leaves = frozenset([l.name for l in tree])
    nodes = [n for n in tree.iter_descendants()]
    descendants = [frozenset(l.name for l in n) for n in nodes]
    mirrors = [leaves-d for d in descendants]

    support = {bp:0 for bp in descendants}
    mirror_support = {bp:0 for bp in mirrors}
    for bst in bs_trees:
        for n in bst.iter_descendants():
            bp = frozenset(l.name for l in n)
            try:
                support[bp] += 1
            except KeyError:
                try:
                    mirror_support[bp] += 1
                except KeyError:
                    pass

    for n in nodes:
        n.support = 0
    for bp, s in support.items():
        idx = descendants.index(bp)
        nodes[idx].support += s
    for bp, s in mirror_support.items():
        idx = mirrors.index(bp)
        nodes[idx].support += s

    support_tree = tree_file.rsplit('.', 1)[0] + '.{}bootstraps.nwk'.format(len(bs_tree_files))
    tree.write(outfile=support_tree)
        
    return support_tree


