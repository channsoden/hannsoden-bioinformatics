#!/usr/bin/env python
# Standard modules
import os, sys, random, argparse, time, itertools, pickle
from argparse import Namespace

# Nonstandard modules
import numpy as np
import pandas as pd
from ete3 import Tree

# My modules
from processing_tools import mapPool
from fasta_tools import get_scaffold_lengths
from fasta_counter import N50
from WholeGenomePhylogeny import main as wgp_main
from WGP2_multiple_alignment import shorten_name
from WGP6_cleanup import cleanup as wgp_cleanup

def parse_args():
    parser = argparse.ArgumentParser(description='Generates a whole genome maximum likelihood phylogeny from a list of fasta genomes.')
    parser.add_argument('-b', '--branch_length', type=float, default=0.333,
                        help=('maximum tolerable branch length in each tree '+
                              '(highly divergent trees tend to have lower quality alignments, and fewer comparable sites)'))
    parser.add_argument('-t', '--tree_size', type=int, default=20, help='maximum number of genomes to compare at a time')
    parser.add_argument('-j', '--jobs', type=int, default=10, help='number of jobs to submit simultaneously')
    parser.add_argument('-f', '--max_effort', type=int, default=20, help='maximum number of passes to make')
    parser.add_argument('-e', '--seed', type=int, default=time.time(), help='set pseudo-random seed for reproducable results')
    parser.add_argument('-o', '--output', type=str, default='genome_clusters', help='basename for output file')
    parser.add_argument('-d', '--dataset', type=str, default='', help='name of a dictionary of fasta genomes in genome_datasets.py')
    parser.add_argument('genomes', type=str, nargs='*', help='list of fasta formatted genomes')
    parser.add_argument('-r', '--restart', dest='restart', action='store_true')
    parser.set_defaults(restart=False)
    args = parser.parse_args()

    if args.restart:
        args = pickle.load(open('genome_clustering.args', 'rb'))
        args.restart = True
        return args
    
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

    pickle.dump(args, open('genome_clustering.args', 'wb'))
    return args

def main():
    global log
    log = sys.stdout
    
    global qualities
    if args.restart:
        qualities = pickle.load(open(args.output+'.N50', 'rb'))
    else:
        qualities = {gen:N50(get_scaffold_lengths(gen).values()) for gen in args.genomes}
        pickle.dump(qualities, open(args.output+'.N50', 'wb'))
    global genome_names
    genome_names = {shorten_name(gen):gen for gen in args.genomes}

    if args.restart:
        clusters = pickle.load(open(args.output+'.clusters', 'rb'))
        completed = pickle.load(open(args.output+'.complete', 'rb'))
    else:
        clusters = [[gen] for gen in args.genomes]
        completed = []
        pickle.dump(clusters, open(args.output+'.clusters', 'wb'))
        pickle.dump(completed, open(args.output+'.complete', 'wb'))
        
    num_genomes = len(args.genomes)
    global comparisons
    if args.restart:
        comparisons = pd.read_csv(args.output+'_progress.tsv', sep='\t', index_col=0)
    else:
        comparisons = pd.DataFrame(np.identity(num_genomes), index=args.genomes, columns=args.genomes)
        comparisons.to_csv(args.output+'_progress.tsv', sep='\t')

    global effort
    effort = 1
    
    while comparisons.sum().sum() != num_genomes **2 and effort <= args.max_effort:
        log.write( 'Began clustering pass #{}\n'.format(effort) )
        groups = random_grouping(clusters)
        jobs = [(compare, (grp, i)) for i, grp in enumerate(groups)]
        results = mapPool(args.jobs, jobs, daemonic=True)
        new_clusters = [clust for clust in results]
        mergers = len(clusters) - len(new_clusters)
        log.write( '{} clusters merged\n'.format(mergers) )
        
        completed += [clust for clust in new_clusters if comparisons[clust[0]].sum() == num_genomes]
        clusters = [clust for clust in new_clusters if clust not in completed]
        effort += 1

        pickle.dump(clusters, open(args.output+'.clusters', 'wb'))
        pickle.dump(completed, open(args.output+'.complete', 'wb'))
        comparisons.to_csv(args.output+'_progress.tsv', sep='\t')

    log.write( '\nFinished clustering\n' )
    total_comparisons = (comparisons.sum().sum() - num_genomes) / 2
    max_comparisons = (num_genomes **2 - num_genomes) / 2
    completion = float(total_comparisons) / max_comparisons * 100.
    log.write( '{}/{} ({}%) genome comparisons made\n'.format(total_comparisons, max_comparisons, completion) )
    log.write( 'clusters = ' )
    log.write( str(clusters)+'\n' )
    log.write( 'completed = ' )
    log.write( str(completed)+'\n' )
                   
    return clusters, completed
        
def best(genomes):
    return max(genomes, key = lambda gen: qualities[gen])

def random_grouping(clusters):
    clusters = clusters[:] # make a copy that we can destroy
    groups = []
    while clusters:
        seed = clusters.pop()
        uncompared = [clust for clust in clusters if not comparisons[clust[0]][seed[0]]]
        random.shuffle(uncompared)
        group = [seed] + uncompared[:args.tree_size-1]
        clusters = [clust for clust in clusters if clust not in group]
        groups.append(group)
    return groups

def compare(group, output_suffix):
    if len(group) < 4:
        # Group is not large enough for meaningful phylogeny
        return group

    short_group = [[shorten_name(name) for name in clust] for clust in group]
    log.write( 'comparing {}\n'.format(short_group) )
    
    best_in_clusters = [best(clust) for clust in group]
    reference = best(best_in_clusters)
    tree_file = wgp(best_in_clusters, reference, output_suffix)
    if not tree_file:
        # Making the tree failed. Abandon this comparison.
        log.write( 'phylogeny of {} failed\n'.format(short_group) )
        try:
            wgp_cleanup(None)
        except:
            pass
        return group
    tree = Tree(tree_file)

    # If the tree has very long branches, than I should split the tree up and compare the members of the subtrees
    # separately. Very long divergences probably reduce the quality and size of the alignment, so strains
    # separated by long branches should be compared separately to utilize more information.
    sub_trees = split_tree(tree)
    if len(sub_trees) > 1:
        # The tree had one or more branches above the threshold set in args.branch_length.
        super_groups = [[genome_names[leaf.name] for leaf in sub_tree.iter_leaves()]
                        for sub_tree in sub_trees]
        # Make sure I never try to compare members of the super groups to each other
        record_comparisons(super_groups)
        # Return the original clusters
        log.write( 'group {} split into {} super groups due to long branch lengths\n'.format(short_group, len(super_groups)) )
        return group
    
    record_comparisons(group)

    threshold = branch_length_threshold(tree)
    if threshold:
        super_clusters = cluster_nodes(tree, threshold)
        log.write( 'group {} clustered into {} super clusters\n'.format(short_group, len(super_clusters)) )
    else:
        # The branches are all very short,
        # the genomes are all very similar to each other.
        # All members of the group are most likely a single cluster
        super_clusters = [group]
        log.write( 'grouping all members of {} due to very short tree\n'.format(short_group) )

    clusters = [merge(super_cluster) for super_cluster in super_clusters]    
    return clusters

def wgp(genomes, reference, output_suffix):
    cwd = os.getcwd()
    output_name = 'phylogeny{}.{}'.format(effort, output_suffix)
    wgp_args = Namespace(reference = reference,
                         step = 1,
                         seed = args.seed,
                         output = output_name,
                         dataset = '',
                         genomes = genomes)
    try: 
        tree_file = wgp_main(wgp_args)
    except IndexError:
        # Probably, these genomes were too divergent to align well, and no alignment was made.
        # Just try again with a different group.
        tree_file = 0
    os.chdir(cwd) # in case wgp_main bugged out and left us in a different directory
    return tree_file

def split_tree(tree):
    # Split the tree into subtrees at long branches.
    long_branches = [n for n in tree.traverse() if n.dist > args.branch_length and not n.is_root()]
    if long_branches:
        breakpoint = long_branches[0]
        treeA = tree
        leavesA = [l for l in breakpoint]
        namesA = [l.name for l in leavesA]
        treeB = tree.copy()
        leavesB = [l for l in treeB if l.name not in namesA]
        if len(leavesA) > 1:
            treeA.prune(leavesA)
        else:
            # breaks when there is only 1 leaf
            # this odd copy construction is necessary to make the single leaf tree a root node
            treeA = Tree(leavesA[0].write())
        if len(leavesB) > 1:
            treeB.prune(leavesB)
        else:
            treeB = Tree(leavesB[0].write())
        return split_tree(treeA) + split_tree(treeB)
    else:
        return [tree]

def record_comparisons(group):
    for i, clusterA in enumerate(group):
        for clusterB in group[i+1:]:
            for genA, genB in itertools.product(clusterA, clusterB):
                comparisons[genA][genB] = 1
                comparisons[genB][genA] = 1

def branch_length_threshold(tree, min_tree_length = 0.05):
    # Determine outlier threshold
    # Threshold, x, is value where
    # x / Q1 = Q1 / Q3
    # Using weighted quartiles
    branch_lengths = [n.dist for n in tree.traverse()]
    branch_lengths.sort()
    total = sum(branch_lengths)
    if total < min_tree_length:
        return 0
    
    divisions = [total / 4., total / 2., total * 3. / 4.]

    Q = []
    cumul = 0
    for l in branch_lengths:
        cumul += l
        try:
            if cumul > divisions[len(Q)]:
                Q.append(l)
        except IndexError:
            break

    #return Q[1] / Q[2] * Q[0] # linient
    return Q[0] ** 2 / Q[2] # strict 

def cluster_nodes(tree, threshold):
    # Find monophyletic clades where all branch lengths are below the threshold
    leaves = [l for l in tree]
    clusters = []
    while leaves:
        last_node = leaves.pop()
        node = last_node.up
        while max([n.dist for n in node.iter_descendants()]) < threshold:
            last_node = node
            node = node.up
        cluster = [l for l in last_node]
        clusters.append([l.name for l in cluster])
        leaves = [l for l in leaves if l not in cluster]
    return clusters

def merge(super_cluster):
    return [gen for clust in super_cluster for gen in clust]

if __name__ == '__main__':
    global args
    args = parse_args()
    random.seed(args.seed)
    clusters, completed = main()


