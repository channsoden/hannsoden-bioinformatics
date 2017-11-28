#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import sys

# Nonstandard modules
from ete3 import Tree

# My modules

"""
def sub_tree(main, tree, args):
    tree = Tree(tree)
    dist, clade1, clade2 = split_tree(tree)

    args1 = reargument(clade1, args)
    args2 = reargument(clade2, args)

    clade1 = main(args1)
    clade2 = main(args2)
    
    return dist, clade1, clade2

def reargument(clade, args):
    if len(clade) > 4:
        # Continue iterating
        pass
    elif len(leaves) > 3:
        # do 1 more iteration
        args.iterate = False
    else:
        # finish
        pass
    return args

    

    sys.exit()
    return args
"""

def split_tree(tree):
    leaves = tree.get_leaves()
    max_dist = 0
    most_distant_pair = ('','')
    for idx, i in enumerate(leaves[:-1]):
        for j in leaves[idx+1:]:
            i2j = tree.get_distance(i, j)
            if i2j > max_dist:       
                max_dist = i2j
                most_distant_pair = (i, j) 

    split_node = most_distant_pair[0].get_common_ancestor(most_distant_pair[1])
    subtrees = split_node.children

    if split_node.is_root():
        split_trees = subtrees
        #split_trees = tuple([subtree.detach() for subtree in subtrees])
    else:
        largest_subtree = max(subtrees, key=lambda subtree: max([subtree.get_distance(l) for l in subtree.iter_leaves()]))
        new_tree = largest_subtree.detach()
        split_trees = (new_tree, tree)

    return split_trees

def break_tree(tree):
    if len(tree.get_leaves()) > 3:
        left, right = split_tree(tree)
        left = break_tree(left)
        right = break_tree(right)
        return '({}, {})'.format(left, right)
    else:
        print tree.write()
        return tree.write().strip(';')
    

tree = Tree(sys.argv[1])
print tree.write()
print break_tree(tree)
