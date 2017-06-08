#!/usr/bin/python

# Load modules.
import math, numpy as np, networkx as nx
import sys, argparse
from collections import defaultdict

from common import combined_similarity_matrix
from hhio import load_index_gene, load_weighted_edge_list, progress

# Parse arguments.
def get_parser():
    description = 'Compute test statistic on hierarchical clusterings.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-elf', '--edge_list_file', type=str, required=True, help='Hierarchy edge list filename')
    parser.add_argument('-igf', '--index_gene_file', type=str, required=True, help='Hierarchy index-gene filename')
    parser.add_argument('-cc', '--cut_criterion', type=str, required=True, choices=['height', 'statistic'])
    parser.add_argument('-ct', '--cut_threshold', type=float, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    return parser

# Define functions.
def statistic(sizes):
    return max(sizes)

def find_cut_height(preordered_T, index_to_gene, threshold, reverse=True):
    '''
    Find clusters for specific height in dendrogram.
    '''
    # Associate nodes of dendrogram with indices.
    index_to_node = defaultdict(set)
    for index, gene in index_to_gene.items():
        index_to_node[index] = frozenset([gene])

    # Initialize clusters with leaf nodes.
    clusters = set(index_to_node.values())

    # Update clusters while ascending dendrogram.
    T = sorted(preordered_T, key=lambda x: x[2], reverse=reverse)
    m = len(T)

    for k, edge in enumerate(T):
        source, target, height = edge
        if (not reverse and height>threshold) or (reverse and height<threshold):
            break

        a = index_to_node[source]
        b = index_to_node[target]
        c = frozenset(set(a) | set(b))

        clusters.discard(a)
        clusters.discard(b)
        clusters.add(c)

        del index_to_node[source]
        index_to_node[target] = c

    return clusters

def find_cut_statistic(preordered_T, index_to_gene, threshold, reverse=True):
    '''
    Find clusters for specific number of nodes in dendrogram.
    '''
    # Associate nodes of dendrogram with indices.
    index_to_node = defaultdict(set)
    for index, gene in index_to_gene.items():
        index_to_node[index] = frozenset([gene])

    # Initialize clusters with leaf nodes.
    clusters = set(index_to_node.values())

    # Update clusters while ascending dendrogram.
    T = sorted(preordered_T, key=lambda x: x[2], reverse=reverse)
    m = len(T)

    for k, edge in enumerate(T):
        source, target, height = edge

        a = index_to_node[source]
        b = index_to_node[target]
        c = frozenset(set(a) | set(b))

        clusters.discard(a)
        clusters.discard(b)
        clusters.add(c)

        del index_to_node[source]
        index_to_node[target] = c

        if statistic(map(len, clusters))>=threshold:
            break

    return clusters

# Run script.
def run(args):
    # Load data.
    if args.verbose:
        progress('Loading data...')

    T = load_weighted_edge_list(args.edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(args.index_gene_file)

    # Load data.
    if args.verbose:
        progress('Processing data...')

    if args.cut_criterion=='height':
        clusters = find_cut_height(T, index_to_gene, args.cut_threshold)
    elif args.cut_criterion=='statistic':
        clusters = find_cut_statistic(T, index_to_gene, args.cut_threshold)

    # Save data.
    if args.verbose:
        progress('Saving data...')

    sorted_clusters = sorted(sorted(map(sorted, clusters)), key=len, reverse=True)
    output_string = '\n'.join('\t'.join(cluster) for cluster in sorted_clusters)

    with open(args.output_file, 'w') as f:
        f.write(output_string)

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
