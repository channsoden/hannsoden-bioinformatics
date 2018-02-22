#!/usr/bin/env python
"""
Implementations of algorithms described in:
Vandal, Conder, & Gentleman (2009) "Minimal Covers of Maximal Cliques for Interval Graphs" Ars Combin 92: 97-129.
"""
# Standard Modules
from itertools import product, chain

# Non-standard Modules
import numpy as np

# My Modules
from processing_tools import mapPool

# Y[:, j] is a vertex in the interval graph
# Y[i, :] is a maximal clique
#                      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
test_data = np.array([[1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], # 0
                      [0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0], # 1
                      [0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0], # 2
                      [0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0], # 3
                      [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0], # 4
                      [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0], # 5
                      [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0], # 6
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1]]) # 7
#                     0  6  9 12 15
reduced_test_data = [[1, 0, 0, 0, 0], # 0
                     [0, 1, 0, 0, 0], # 1
                     [0, 1, 1, 0, 0], # 2
                     [0, 1, 1, 0, 0], # 3
                     [0, 0, 1, 0, 0], # 4
                     [0, 0, 0, 1, 0], # 5
                     [0, 0, 0, 1, 0], # 6
                     [0, 0, 0, 0, 1]] # 7


# Solution includes
[[0, 7, 2, 5],
 [0, 7, 2, 6],
 [0, 7, 3, 5],
 [0, 7, 3, 6]]

#                    0  1  2  3  4  5  6
my_test = np.array([[1, 1, 1, 0, 0, 0, 0], # 0
                    [0, 1, 1, 1, 0, 0, 0], # 1
                    [0, 0, 0, 1, 1, 0, 0], # 2
                    [0, 0, 0, 0, 1, 1, 0], # 3
                    [0, 0, 0, 0, 0, 1, 1]]) # 4

split_test = np.array([[1, 0, 0, 0, 0, 0],
                       [0, 1, 0, 0, 0, 0],
                       [0, 0, 1, 0, 0, 0],
                       [0, 0, 1, 1, 0, 0],
                       [0, 0, 1, 1, 0, 0],
                       [0, 0, 0, 1, 0, 0],
                       [0, 0, 0, 0, 1, 1]])

# Solution is
[[0, 2, 4],
 [0, 1, 3, 4]]

def enumerate_minimal_covers(clique_matrix, threads = 1):
    m, n = clique_matrix.shape
    elements = get_duals(clique_matrix, n)
    elements = reduce_elements(elements, m, n)
    elements = order_minimals(elements)
    islands = split_disconnected(elements)
    parts = [(minimal_covers, island) for island in islands]
    partial_covers, partial_chains = list(zip(*mapPool(threads, parts)))
    covers = merge_disconnected(partial_covers)
    return covers

class vertex(object):
    def __init__(self, j, cliques):
        self.j = j
        self.dual = np.nonzero(cliques)[0]
        self.min = np.min(self.dual)
        self.max = np.max(self.dual)
        self.len = len(self.dual)

    def __repr__(self):
        return 'v.{}'.format(self.j)

class counter(object):
    def __init__(self, start):
        self.count = start
    def tick(self):
        self.count += 1
        if not self.count % 10000:
            print('. . . {} recursive calls made'.format(self.count))
    
def get_duals(Y, n):
    """Y is a clique matrix.
    Rows (axis 1) in Y represent maximal cliques of the interval graph.
    Columns (axis 0) in Y represent subset-minimal vertices of the interval graph.
    Values of Y are either 1, denoting inclusion of the vertex in the clique, and the clique in the dual of the vertex."""
    duals = [vertex(j, Y[:, j]) for j in range(n)]
    return set(duals)
    
def reduce_elements(elements, m, n):
    """Reduce a set of vertex objects to only subset-minimal elements.
    A vertex is subset-minimal if it's dual has no other dual as a proper subset."""
    
    # First, if two elements share an endpoint, keep only the shortest one
    elements_by_min = [('none', None)] * m # str > int, will sort to end (None will sort to front)
    elements_by_max = [('none', None)] * m
    for d in elements:
        elements_by_min[d.min] = min([(d.len, d), elements_by_min[d.min]])
        elements_by_max[d.max] = min([(d.len, d), elements_by_max[d.max]])
    shortest_by_min = set(zip(*elements_by_min)[1])
    shortest_by_max = set(zip(*elements_by_max)[1])
    elements = shortest_by_min & shortest_by_max & elements

    # Keep ABAB patterns, toss A in ABBA patterns
    endpoints = []
    for d in elements:
        endpoints.append((d.min, d))
        endpoints.append((d.max, d))
    endpoints.sort()

    supersets = set([])
    open_elements = []
    for pos, d in endpoints:
        if d in open_elements:
            supersets |= set(open_elements[:open_elements.index(d)])
            open_elements.remove(d)
        else:
            open_elements.append(d)
    elements -= supersets
    
    return elements

def order_minimals(reduced):
    reduced = list(reduced)
    reduced.sort(key= lambda v: v.min)
    return reduced

def minimal_covers(minimals, m, last_clique=[], last_element=[], effort=None):
    """Returns a list of all minimal covers of an interval graph and a corresponding chain for each minimal cover.
    minimals is the set of subset-minimals ordered by their left/right endpoints.
    m is the number of maximal cliques in the interval graph.
    last_clique is the index of the last clique selected by the algorithm.
    last_element is the previous iteration's minimals[0]

    Complexity of this algorithm is O(2**(m-2)). Should be feasible for m <= 30 ish."""

    if not last_clique:
        # Top level
        if m > 30:
            effort = counter(1)
            print('attempting to enumerate minimals covers with m = {}'.format(m))
            
        minimals, simplicial_elements, essential_cliques = remove_simplicial(minimals)

        if not minimals:
            # There is a single minimal cover
            return [essential_cliques], [simplicial_elements]
    else:
        essential_cliques = simplicial_elements = []
        if effort: effort.tick()
        
    if not minimals:
        # Bottom level
        return [last_clique], [last_element] # changed here, included last_clique
    else:
        covers = []
        chains = []
        
        first_right_end = minimals[0] # L1
        # can use first minimal since minimals is both left and right end ordered.
        # argmin(max(y*)) will always be the first minimal.

        # clique selection should be difference of first_right_end and the previous first_right_end
        clique_selection = [i for i in first_right_end.dual if (not last_element or last_element[0].max < i)] # L2 (set to loop through)
        
        # to return only the minimum covers on the test data
        # only keep the mus that cover the maximum number of minimals
        # but I can't prove it will always work

        for clique in clique_selection: # L2
            remaining_minimals = [v for v in minimals if clique < v.min] # R1 (first argument)
            # R1, changed here, pass both last_clique and last_element
            partial_covers, partial_chains = minimal_covers(remaining_minimals, m,
                                                            last_clique=[clique], last_element=[first_right_end],
                                                            effort=effort) 
            for cover in partial_covers: # L3
                covers.append(essential_cliques + last_clique + cover) # changed here, add last_clique rather than last_element.max
            for chain in partial_chains:
                chains.append(simplicial_elements + last_element + chain)

    return covers, chains

def remove_simplicial(elements):
    # Observation 3.6.1
    # Remove simplicial elements, and a priori select the maximal cliques that contain them.
    simplicial = [v for v in elements if v.len == 1]
    essential_cliques = [v.min for v in simplicial]
    reduced = [v for v in elements if v not in simplicial]
    return reduced, simplicial, essential_cliques

def split_disconnected(elements):
    """Observation 3.6.2
    The problem can be broken down into parts, doing each connected component of the interval graph separately, then combining the results.
    The full result is thus the set of all combinations of partial results.

    This function requires elements to be subset-minimals sorted by their left/right endpoints."""
    i = 0
    islands = [[]]
    sizes = []
    for element in elements:
        if element.min <= i:
            islands[-1].append(element)
        else:
            sizes.append(1 + islands[-1][-1].max - islands[-1][0].min)
            islands.append([element])
        i = element.max
    sizes.append(1 + islands[-1][-1].max - islands[-1][0].min)
    return list(zip(islands, sizes))

def merge_disconnected(partial_covers):
    """Combines partial results for each disconnected region of the interval graph into the full set of all covers."""
    grouped_combos = product(*partial_covers)
    combos = [list(chain(*grouped_combo)) for grouped_combo in grouped_combos]
    return combos


if __name__ == '__main__':
    print(enumerate_minimal_covers(test_data))
