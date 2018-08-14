#!/usr/bin/python

import math, numpy as np

####################################################################################################################################
#
# General IO functions
#
####################################################################################################################################

def load_index_gene(filename):
    '''
    Load index-gene associations.
    '''
    index_to_gene = dict()
    gene_to_index = dict()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.strip().split()
                index = int(arrs[0])
                gene = arrs[1]
                index_to_gene[index] = gene
                gene_to_index[gene] = index
    return index_to_gene, gene_to_index

def save_index_gene(filename, index_to_gene):
    '''
    Save index-gene associations.
    '''
    index_gene_list = sorted(index_to_gene.items())
    with open(filename, 'w') as f:
        f.write('\n'.join('{}\t{}'.format(index, gene) for index, gene in index_gene_list))

def load_edge_list(filename, dictionary=None, unweighted=False):
    '''
    Load edge list.
    '''
    edge_list = list()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.strip().split()

                if len(arrs)==2:
                    i = int(arrs[0])
                    j = int(arrs[1])
                    edge_list.append((i, j, 1))
                elif len(arrs)==3:
                    i = int(arrs[0])
                    j = int(arrs[1])
                    try:
                        weight = float(arrs[2])
                        if not math.isinf(weight) and not math.isnan(weight):
                            edge_list.append((i, j, weight))
                        else:
                            raise Warning('{} is not a valid edge weight; edge omitted.'.format(arrs[2]))
                    except ValueError:
                        raise Warning('{} is not a valid edge weight; edge omitted.'.format(arrs[2]))
                elif l.strip():
                    raise Warning('{} is not a valid edge; edge omitted.'.format(l.strip()))

    if dictionary:
        edge_list = [(dictionary[i], dictionary[j], weight) for i, j, weight in edge_list]

    if unweighted:
        edge_list = [(i, j) for i, j, weight in edge_list]

    return edge_list

def save_edge_list(filename, edge_list, dictionary=None):
    '''
    Save edge list.
    '''
    if dictionary:
        edge_list = [(dictionary[edge[0]], dictionary[edge[1]]) + tuple(edge)[2:] for edge in edge_list]

    with open(filename, 'w') as f:
        f.write('\n'.join('\t'.join(map(str, edge)) for edge in edge_list))

def load_gene_score(filename, score_threshold=0.0):
    '''
    Load gene scores.
    '''
    gene_to_score = dict()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.strip().split()
                if len(arrs)==2:
                    gene = arrs[0]
                    try:
                        score = float(arrs[1])
                        if not math.isinf(score) and not math.isnan(score):
                            if score>=score_threshold:
                                gene_to_score[gene] = score
                        else:
                            raise Warning('{} is not a valid gene score; gene score omitted.'.format(arrs[1]))
                    except ValueError:
                        raise Warning('{} is not a valid gene score; gene score omitted.'.format(arrs[1]))
                elif l.strip():
                    raise Warning('{} is not a valid gene score; gene score omitted.'.format(l.strip()))

    return gene_to_score

def save_gene_score(filename, gene_to_score):
    '''
    Save gene scores.
    '''
    gene_score_list = sorted(gene_to_score.items())
    with open(filename, 'w') as f:
        f.write('\n'.join('{}\t{}'.format(gene, score) for gene, score in gene_score_list))

####################################################################################################################################
#
# Matrix IO functions
#
####################################################################################################################################

def load_matrix(filename, matrix_name='PPR'):
    '''
    Load matrix.
    '''
    import h5py

    f = h5py.File(filename, 'r')
    if matrix_name in f:
        A = np.asarray(f[matrix_name].value, dtype=np.float32)
    else:
        raise KeyError('Matrix {} is not in {}.'.format(matrix_name, filename))
    f.close()
    return A

def save_matrix(filename, A, matrix_name='PPR'):
    '''
    Save matrix.
    '''
    import h5py

    f = h5py.File(filename, 'a')
    if matrix_name in f:
        del f[matrix_name]
    f[matrix_name] = np.asarray(A, dtype=np.float32)
    f.close()

####################################################################################################################################
#
# Screen IO functions
#
####################################################################################################################################

def progress(message=''):
    '''
    Write status message to screen; overwrites previous message and does not advance line.
    '''
    import sys

    try:
        rewind = progress.rewind
    except AttributeError:
        rewind = 0

    sys.stdout.write('\r'+' '*rewind)
    sys.stdout.flush()
    sys.stdout.write('\r'+str(message))
    sys.stdout.flush()
    progress.rewind = max(len(str(message).expandtabs()), rewind)
