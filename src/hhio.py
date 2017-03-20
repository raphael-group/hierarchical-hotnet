#!/usr/bin/python

from common import convert_edge_list, convert_weighted_edge_list

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

def load_edge_list(filename, dictionary=None):
    '''
    Load edge list.
    '''
    edge_list = list()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.strip().split()
                i = int(arrs[0])
                j = int(arrs[1])
                edge_list.append((i, j))

    if dictionary:
        edge_list = convert_edge_list(edge_list, dictionary)

    return edge_list

def save_edge_list(filename, edge_list, dictionary=None):
    '''
    Save edge list.
    '''
    if dictionary:
        edge_list = convert_edge_list(edge_list, dictionary)

    with open(filename, 'w') as f:
        f.write('\n'.join('{}\t{}'.format(i, j) for i, j in edge_list))

def load_gene_score(filename, score_threshold=0.0):
    '''
    Load gene scores.
    '''
    gene_to_score = dict()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.strip().split()
                gene = arrs[0]
                score = float(arrs[1])
                if score>=score_threshold:
                    gene_to_score[gene] = score
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
        A = f[matrix_name].value
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
    f[matrix_name] = A
    f.close()

####################################################################################################################################
#
# Hierarchy IO functions
#
####################################################################################################################################

def load_weighted_edge_list(filename, dictionary=None):
    '''
    Load weighted edge list from a TSV file in (source, target, weight) format.
    '''
    edge_list = list()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.rstrip('\n').split('\t')
                source = int(arrs[0])
                target = int(arrs[1])
                weight = float(arrs[2])
                edge_list.append((source, target, weight))

    if dictionary:
        edge_list = convert_weighted_edge_list(edge_list, dictionary)

    return edge_list

def save_weighted_edge_list(filename, edge_list, dictionary=None):
    '''
    Save weighted edge list as a TSV file in (source, target, weight) format.
    '''
    if dictionary:
        edge_list = convert_weighted_edge_list(edge_list, dictionary)

    with open(filename, 'w') as f:
        f.write('\n'.join('\t'.join(map(str, edge)) for edge in edge_list))

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
