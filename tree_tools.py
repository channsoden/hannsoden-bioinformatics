#!/usr/bin/env python
# Standard modules
import sys

# Nonstandard modules
from ete3 import Tree

def root(tree, outgroup):
    if len(outgroup) == 1:
        tree.set_outgroup(outgroup[0])
    else:
        ancestor = tree.get_common_ancestor(*outgroup)
        tree.set_outgroup(ancestor)
    
def clone_correct(tree, genomes):
    pass

def change_names(tree, name_dict):
    for leaf in tree.iter_leaves():
        leaf.name = name_dict[leaf.name]

if __name__ == '__main__':
    from genome_datasets import summer2016_neurospora
    
    tree_file = sys.argv[1]
    tree = Tree(tree_file)

    change_names(tree, {v:k for k, v in summer2016_neurospora.items()})
    simplify(tree, tree_file+'.simple')
