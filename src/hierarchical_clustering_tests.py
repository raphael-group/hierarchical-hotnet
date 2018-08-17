#!/usr/bin/python

import math, numpy as np

from hierarchical_clustering import *
from hhio import progress

####################################################################################################################################
#
# Test functions
#
####################################################################################################################################

def naive_HD(A, reverse=False):
    '''
    Compute the hierarchical decomposition of a graph into its strongly connected components using a naive approach.
    '''
    import networkx as nx

    # Collect nodes and edges.
    num_nodes = np.shape(A)[0]
    nodes = list(range(num_nodes))
    generator = ((i, j, A[i, j]) for i in nodes for j in nodes if i!=j and A[i, j]!=0)
    edges = sorted(generator, reverse=reverse, key=lambda edge: edge[2])
    num_edges = len(edges)

    # Initialize hierarchy.
    num_components = num_nodes
    j = num_nodes
    ancestors = nodes
    T = list()

    # Build hierarchy by adding edges to graph ordered by weight.
    G = nx.DiGraph()
    G.add_nodes_from(nodes)

    for k in range(num_edges):
        G.add_edge(edges[k][0], edges[k][1], weight=edges[k][2])
        if k==num_edges-1 or edges[k][2]!=edges[k+1][2]:
            components = list(nx.strongly_connected_components(G))
            if len(components)!=num_components:
                num_components=len(components)
                for component in components:
                    component_ancestors = set(ancestors[i] for i in component)
                    if len(component_ancestors)>1:
                        for i in component:
                            ancestors[i] = j
                        for i in component_ancestors:
                            T.append((i, j, edges[k][2]))
                        j += 1
            if len(components)==1:
                break

    T = [(u+1, v+1, w) for u, v, w in T]

    return T

def edges_to_adjacency_matrix(edges):
    '''
    Create adjacency matrix from edge list.
    '''
    nodes = sorted(set(edge[0] for edge in edges) | set(edge[1] for edge in edges))
    num_nodes = len(nodes)
    node_to_index = dict((node, i) for i, node in enumerate(nodes))

    weight_type = set(type(edge[2]) for edge in edges if len(edge)>2)
    if weight_type==set([int]):
        dtype = np.int
    else:
        dtype = np.float64

    A = np.zeros((num_nodes, num_nodes), dtype=dtype)
    for edge in edges:
        if len(edge)==2:
            u, v = edge
            w = 1
        elif len(edge)==3:
            u, v, w = edge

        A[node_to_index[u], node_to_index[v]] = w

    return nodes, A

def tarjan_1983_example():
    """
    Generate the graph from the example in Figure 1 of Tarjan (1983).
    """
    E = [('a', 'b', 10), ('b', 'a', 12), ('b', 'c', 30), ('d', 'c',  6), ('d', 'e', 16),
         ('e', 'd', 13), ('e', 'f',  8), ('f', 'a', 26), ('a', 'g', 15), ('g', 'b', 35),
         ('c', 'g', 45), ('g', 'c', 22), ('d', 'g', 14), ('g', 'e', 50), ('f', 'g', 20)]

    V, A = edges_to_adjacency_matrix(E)

    return V, A, E

def random_adjacency_matrix(n,seed=np.random.randint(0,4294967295),sparsity=0.0,nonuniqueness=0.0):

    np.random.seed(seed=seed)

    sparsity = max(0,min(sparsity,1))
    nonuniqueness = max(0,min(nonuniqueness,1))
    unique_elements = int((1-sparsity)*(1-nonuniqueness)*n**2)
    repeated_elements = int((1-sparsity)*nonuniqueness*n**2)

    # Should the graph be complete with unique edge weights?

    if sparsity==0 and nonuniqueness==0:

        # If so, generate a random matrix and remove the diagonal entries.

        A = np.random.rand(n,n)
        np.fill_diagonal(A,0)

    else:

        # If not, generate a vector with the expected number of unique
        # entries, duplicate some of those entries according to a Binomial
        # distribution according to the nonuniqueness parameter, and fill
        # in the remaining entries with zeros according to the sparsity
        # parameter.  Permute the entries of the vector and reshape it into
        # the correct dimensions.

        B = np.random.rand(unique_elements)
        C = np.zeros(n**2-unique_elements,dtype=np.float)

        tally = 0
        while tally<repeated_elements:
            duplications = max(np.random.binomial(repeated_elements-tally,0.5),1)
            C[tally:tally+duplications] = B[np.random.randint(unique_elements)]
            tally += duplications

        A = np.random.permutation(np.concatenate((B,C))).reshape(n,n)
        np.fill_diagonal(A,0)

        # Is the resulting graph strongly connected?  If not, add entries
        # to the adjacency matrix at semi-random to bridge the components.

        components = strongly_connected_components(A)
        m = len(components)
        if m>1:
            for i in range(m):
                for j in range(m):
                        if i!=j:
                            p = np.random.randint(len(components[i]))
                            q = np.random.randint(len(components[j]))
                            A[components[i][p],components[j][q]] = np.random.rand()

    return A

def are_trees_equivalent(S, T, index_to_gene, reverse=False):
    '''
    Check whether two hierarchies are equivalent, i.e., whether they produce the same partitions of
    the leaf nodes at the same heights.
    '''
    height_to_clusters_S = find_height_to_clusters(S, index_to_gene, reverse=reverse)
    height_to_clusters_T = find_height_to_clusters(T, index_to_gene, reverse=reverse)
    return height_to_clusters_S==height_to_clusters_T

def show_hierarchy(T, reverse=False):
    y = list()
    for source, target, height in sorted(T, key=lambda x: (x[2], x[1], x[0] if not reverse else -x[0])):
        y.append('\t{}\t{}\t{}'.format(source, target, height))
    return '\n'.join(y)

####################################################################################################################################
#
# Tests
#
####################################################################################################################################

def examples():
    import time

    # Example 1
    print('Example 1: Tarjan (1983)')
    reverse = False

    # Generate example from Tarjan (1983).
    V, A, E = tarjan_1983_example()
    index_to_gene = dict((i+1, v) for i, v in enumerate(V))

    # Compute hierarchy using implementation of naive algorithm.
    S = naive_HD(A, reverse=reverse)

    print('\n- Output from naive algorithm')
    print(show_hierarchy(S, reverse=reverse))

    # Compute hierarchy using implementation of Tarjan's algorithm.
    T = tarjan_HD(A, reverse=reverse)

    print('\n- Output from Tarjan\'s algorithm')
    print(show_hierarchy(T, reverse=reverse))

    # Check results.
    if are_trees_equivalent(S, S, index_to_gene, reverse):
        print('\n- The results from the naive algorithm are self-consistent.')
    else:
        print('- The results from the naive algorithm are NOT self-consistent.')
    if are_trees_equivalent(T, T, index_to_gene, reverse):
        print('- The results from Tarjan\'s algorithm are self-consistent.')
    else:
        print('- The results from Tarjan\'s algorithm are NOT self-consistent.')
    if are_trees_equivalent(S, T, index_to_gene, reverse):
        print('- The results from the naive algorithm and Tarjan\'s algorithm are consistent.')
    else:
        print('- The results from the naive algorithm and Tarjan\'s algorithm are NOT consistent.')

    # Example 2
    print('\nExample 2: Tarjan (1983), edges added in reverse')
    reverse = True

    # Generate example from Tarjan (1983).
    V, A, E = tarjan_1983_example()
    index_to_gene = dict((i+1, v) for i, v in enumerate(V))

    # Compute hierarchy using implementation of naive algorithm.
    S = naive_HD(A, reverse=reverse)

    print('\n- Output from naive algorithm')
    print(show_hierarchy(S, reverse=reverse))

    # Compute hierarchy using implementation of Tarjan's algorithm.
    T = tarjan_HD(A, reverse=reverse)

    print('\n- Output from Tarjan\'s algorithm')
    print(show_hierarchy(T, reverse=reverse))

    # Check results.
    if are_trees_equivalent(S, S, index_to_gene, reverse):
        print('\n- The results from the naive algorithm are self-consistent.')
    else:
        print('- The results from the naive algorithm are NOT self-consistent.')
    if are_trees_equivalent(T, T, index_to_gene, reverse):
        print('- The results from Tarjan\'s algorithm are self-consistent.')
    else:
        print('- The results from Tarjan\'s algorithm are NOT self-consistent.')
    if are_trees_equivalent(S, T, index_to_gene, reverse):
        print('- The results from the naive algorithm and Tarjan\'s algorithm are consistent.')
    else:
        print('- The results from the naive algorithm and Tarjan\'s algorithm are NOT consistent.')

    # Example 3
    print('\nTest 3: Various random matrices')
    at = time.time()
    bt = time.time()
    ct = time.time()

    naive_elapsed_time = 0.0
    tarjan_elapsed_time = 0.0

    num_trials = 0
    num_failures = 0
    for reverse in [True, False]:
        for sparsity in [0.0, 0.1, 0.5]:
            for nonuniqueness in [0.0, 0.1, 0.5]:
                for n in [2, 3, 5, 8, 13, 21, 34, 55, 89]:
                    for seed in range(1, 4):
                        progress('- Running trials: n: {}, sparsity: {}, nonuniqueness: {}, reverse: {}, seed: {}...'.format(n, sparsity, nonuniqueness, reverse, seed))

                        A = random_adjacency_matrix(n, seed=seed, sparsity=sparsity, nonuniqueness=nonuniqueness)
                        index_to_gene = dict((i+1, v) for i, v in enumerate(range(n)))

                        at = time.time()
                        S = naive_HD(A, reverse=reverse)
                        bt = time.time()
                        T = tarjan_HD(A, reverse=reverse)
                        ct = time.time()

                        naive_elapsed_time += bt-at
                        tarjan_elapsed_time += ct-bt

                        num_trials += 1
                        if not are_trees_equivalent(S, T, index_to_gene, reverse):
                            num_failures += 1
                            print('Failure.')
                            print('Naive algorithm:')
                            print(show_hierarchy(S))
                            print('Tarjan\'s algorithm:')
                            print(show_hierarchy(T))

    progress()

    print('\n- Number of trials: {}'.format(num_trials))
    print('- Number of failures: {}'.format(num_failures))


    print('- Total elapsed time for naive algorithm: {:.1f} seconds'.format(naive_elapsed_time))
    print('- Total elapsed time for Tarjan\'s algorithm: {:.1f} seconds'.format(tarjan_elapsed_time))

    print('\nTest 4: Various random matrices')
    print('')
    print('- Elapsed time in seconds for Tarjan\'s algorithm:')
    print('\t{}\t{}\t{}\t{}\t{}'.format('Nodes', 'Forward     ', 'Order', 'Backward    ', 'Order'))

    for i in range(5):
        n = 100*2**i

        A = random_adjacency_matrix(n, seed=i, sparsity=0.1, nonuniqueness=0.1)
        at = time.time()
        S = tarjan_HD(A, reverse=False)
        bt = time.time()
        T = tarjan_HD(A, reverse=True)
        ct = time.time()

        if i>0:
            print('\t{}\t{:.2e}\t{:.1f}\t{:.2e}\t{:.1f}'.format(n, bt-at, math.log((bt-at)/previous_bt_at)/math.log(2), ct-bt, math.log((ct-bt)/previous_ct_bt)/math.log(2)))
        else:
            print('\t{}\t{:.2e}\t{:.1f}\t{:.2e}\t{:.1f}'.format(n, bt-at, float('nan'), ct-bt, float('nan')))

        previous_bt_at = bt-at
        previous_ct_bt = ct-bt

if __name__ == '__main__':
    examples()
