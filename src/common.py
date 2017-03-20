#!/usr/bin/python

import math, numpy as np

####################################################################################################################################
#
# Method functions
#
####################################################################################################################################

def degree_sequence(A):
    '''
    Find the degree sequence for an adjacency matrix.
    '''
    return np.sum(A, axis=0)

def laplacian_matrix(A):
    '''
    Find the Laplacian matrix.
    '''
    return np.diag(degree_sequence(A))-A

def walk_matrix(A):
    '''
    Find the walk matrix for the random walk.
    '''
    d = degree_sequence(A)
    d[np.where(d<=0)] = 1
    return np.asarray(A, dtype=np.float64)/np.asarray(d, dtype=np.float64)

def hotnet_similarity_matrix(A, t):
    '''
    Perform the diffusion process in HotNet.
    '''
    from scipy.linalg import eigh
    D, V = eigh(-t*laplacian_matrix(A))
    return np.dot(np.exp(D)*V, sp.linalg.inv(V))

def hotnet2_similarity_matrix(A, beta):
    '''
    Perform the random walk with restart process in HotNet2.
    '''
    from scipy.linalg import inv
    return beta*inv(np.eye(*np.shape(A))-(1-beta)*walk_matrix(A))

def hotnet2_similarity_matrix_iterative(A, beta, steps=1):
    '''
    Perform the random walk with restart process in HotNet2 iteratively.
    '''
    W = walk_matrix(A)
    P = np.eye(*np.shape(A))
    F = np.eye(*np.shape(A))
    for _ in range(steps):
        P[:, :] = (1-beta)*np.dot(W, P)+beta*F
    return P

def hh_similarity_matrix(A, beta):
    '''
    Construct the Hierarchical HotNet similarity matrix.
    '''
    return hotnet2_similarity_matrix(A, beta)

def combined_similarity_matrix(P, gene_to_index, gene_to_score):
    '''
    Find the combined similarity submatrix from the topological similarity matrix and the gene
    scores.
    '''
    m, n = np.shape(P)
    min_index = min(gene_to_index.values())
    max_index = max(gene_to_index.values())

    # We do not currently support nodes without edges, which seems to be the most likely way for
    # someone to encounter this error.
    if m!=n:
        raise IndexError('The similarity matrix is not square.')
    if max_index-min_index+1!=n:
        raise IndexError('The similarity matrix and the index-gene associations have different dimensions.')

    topology_genes = set(gene_to_index)
    score_genes = set(gene_to_score)
    common_genes = sorted(topology_genes & score_genes)
    common_indices = [gene_to_index[gene]-min_index for gene in common_genes]

    f = np.array([gene_to_score[gene] for gene in common_genes])
    S = P[np.ix_(common_indices, common_indices)]*f
    common_index_to_gene = dict((i+1, gene) for i, gene in enumerate(common_genes))
    common_gene_to_index = dict((gene, i) for i, gene in common_index_to_gene.items())

    return S, common_index_to_gene, common_gene_to_index

####################################################################################################################################
#
# Data structure functions
#
####################################################################################################################################

def convert_edge_list_to_adjacency_matrix(edge_list):
    '''
    Convert an edge list to an adjacency matrix for an undirected graph.
    '''
    k = min(min(edge) for edge in edge_list)
    l = max(max(edge) for edge in edge_list)

    A = np.zeros((l-k+1, l-k+1), dtype=np.int)
    for i, j in edge_list:
        A[i-k, j-k] = A[j-k, i-k] = 1

    return A

def convert_adjacency_matrix_to_edge_list(A):
    '''
    Convert an adjacency matrix for an undirected graph to an edge list.
    '''
    m, n = np.shape(A)
    return ((i, j) for i in range(m) for j in range(i, n) if A[i, j]!=0)

def convert_edge_list(edge_list, dictionary):
    '''
    Convert an edge list with indices/names to an edge list with names/indices.
    '''
    return [(dictionary[i], dictionary[j]) for i, j in edge_list]

def convert_weighted_edge_list(edge_list, dictionary):
    '''
    Convert an edge list with indices/names to an edge list with names/indices.
    '''
    return [(dictionary[i], dictionary[j], w) for i, j, w in edge_list]

####################################################################################################################################
#
# Statistics functions
#
####################################################################################################################################

def multiple_hypothesis_correction(p_values_, method='BH'):
    """
    Compute a multiple-hypothesis correction using one of multiple
    methods: Bonferroni (bonferroni), Holm-Bonferroni (holm),
    Benjamini-Hochberg (BH), or Benjamini-Hochberg-Yekutieli (BY).

    Example (from Benjamini and Hochberg (1995)):
    In [1]: p_values = [0.0201, 0.7590, 0.0344, 0.0278, 0.0095,
                        0.0298, 0.6719, 1.0000, 0.0001, 0.4262,
                        0.0004, 0.6528, 0.3240, 0.0459, 0.0019]
    In [2]: q_values = multiple_hypothesis_correction(p_values, method='BH')
    In [3]: q_values
    Out[3]: array([ 0.0603    ,  0.81321429,  0.0645    ,  0.06385714,  0.035625  ,
                    0.06385714,  0.77526923,  1.        ,  0.0015    ,  0.58118182,
                    0.003     ,  0.77526923,  0.486     ,  0.0765    ,  0.0095    ])
    """
    if method not in ['bonferroni', 'holm', 'BH', 'BY']:
        raise NotImplementedError('{} method not implemented'.format(method))

    valid_indices = [i for i, p_value in enumerate(p_values_) if 0.0<=p_value<=1.0]
    invalid_indices = [i for i, p_value in enumerate(p_values_) if not 0.0<=p_value<=1.0]

    p_values = np.asarray(p_values_)[valid_indices]
    n = len(p_values)

    if method=='bonferroni':
        q_values = np.minimum(n*p_values, 1)

    elif method=='holm':
        sorted_indices = np.argsort(p_values)
        sorted_p_values = p_values[sorted_indices]

        sorted_q_values = np.zeros(n)
        sorted_q_values[0] = min(sorted_p_values[0], 1.0)
        for i in range(1, n-1):
            sorted_q_values[i] = min(max(float(n-i)*sorted_p_values[i], sorted_q_values[i-1]), 1.0)

        q_values = np.zeros(n)
        q_values[sorted_indices] = sorted_q_values

    elif method=='BH':
        sorted_indices = np.argsort(p_values)
        sorted_p_values = p_values[sorted_indices]

        sorted_q_values = np.zeros(n)
        sorted_q_values[n-1] = min(sorted_p_values[n-1], 1.0)
        for i in reversed(range(n-1)):
            sorted_q_values[i] = min(float(n)/float(i+1)*sorted_p_values[i], sorted_q_values[i+1])

        q_values = np.zeros(n)
        q_values[sorted_indices] = sorted_q_values

    elif method=='BY':
        sorted_indices = np.argsort(p_values)
        sorted_p_values = p_values[sorted_indices]

        c = np.sum(1.0/np.arange(1, n+1, dtype=np.float64))
        sorted_q_values = np.zeros(n)
        sorted_q_values[n-1] = min(sorted_p_values[m-1], 1.0)
        for i in reversed(range(n-1)):
            sorted_q_values[i] = min(float(n)/float(i+1)*sorted_p_values[i], sorted_q_values[i+1])

        q_values = np.zeros(n)
        q_values[sorted_indices] = sorted_q_values

    q_values_ = np.zeros(len(p_values_))
    q_values_[valid_indices] = q_values
    q_values_[invalid_indices] = float('nan')

    return q_values_
