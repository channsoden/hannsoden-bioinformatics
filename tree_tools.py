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

def tree_order(reference, tree, root=None, include_reference=False):
    # Returns list of leaf names in order of ascending distance away from the reference leaf.
    # reference (str) - leaf name.
    # tree (ete3.TreeNode)
    # root (ete3.TreeNode)
    # include_reference (bool) - if True, reference will be first element in returned list.
    if not root:
        root = tree.get_tree_root()
    order = [reference]
    node = tree&reference
    while True:
        leaves = [l.name for l in node.get_leaves() if not l.name in order]
        if len(leaves) > 2:
            # need to get subtree order
            subtree = [child for child in node.children if not tree&reference in child.get_leaves()][0]
            leaves = tree_order(leaves[0], subtree, root=subtree, include_reference=True)
        order.extend(leaves)
        if node == root:
            break
        node = node.up
    if not include_reference:
        order = order[1:]
    return order

def traverse_between(A, B, include_ancestor=True):
    # Returns list of all nodes between A and B in order
    tree = A.get_tree_root()
    ancestor = tree.get_common_ancestor(A, B)
    
    ascending = []
    node = A
    while node != ancestor:
        ascending.append(node)
        node = node.up
    descending = []
    node = B
    while node != ancestor:
        descending.append(node)
        node = node.up
    descending.reverse()

    if include_ancestor:
        between = ascending + [ancestor] + descending
    else:
        between = ascending + descending
    return between
        
if __name__ == '__main__':
    from genome_datasets import summer2016_neurospora
    
    tree_file = sys.argv[1]
    tree = Tree(tree_file)

    change_names(tree, {v:k for k, v in summer2016_neurospora.items()})
    simplify(tree, tree_file+'.simple')
