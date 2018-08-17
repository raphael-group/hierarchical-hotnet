#!/usr/bin/python

# Load modules.
import numpy as np
import sys, argparse
from collections import defaultdict

from hhio import load_index_gene, load_edge_list, progress

# Parse arguments.
def get_parser():
    description = 'Find cluster sizes across cuts of hierarchy.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-elf', '--edge_list_file', type=str, required=True, help='Hierarchy edge list filenames')
    parser.add_argument('-igf', '--index_gene_file', type=str, required=True, help='Hierarchy index-gene filenames')
    parser.add_argument('-csf', '--cluster_size_file', type=str, required=True, help='Cluster size file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    return parser

# Define functions.
def compute_statistic(sizes):
    return max(sizes)

def compute_statistics(T, index_to_gene, reverse=True):
    heights = list()
    statistics = list()

    # Associate nodes of dendrogram with indices.
    index_to_node = defaultdict(set)
    for index, gene in index_to_gene.items():
        index_to_node[index] = set([gene])

    # Add cluster statistic of leaf nodes.
    height = float('inf') if reverse else 0.0
    sizes = (len(node) for node in index_to_node.values())
    statistic = compute_statistic(sizes)
    heights.append(height)
    statistics.append(statistic)

    # Add cluster statistics for inner nodes.
    T = sorted(T, key=lambda x: x[2], reverse=reverse)
    m = len(T)

    for i, edge in enumerate(T):
        source, target, height = edge

        a = index_to_node[source]
        b = index_to_node[target]
        c = set(a) | set(b)

        del index_to_node[source]
        index_to_node[target] = c

        if i==m-1 or height!=T[i+1][2]:
            sizes = (len(node) for node in index_to_node.values())
            statistic = compute_statistic(sizes)
            heights.append(height)
            statistics.append(statistic)

    # Add cluster statistic for root node.
    height = 0.0 if reverse else float('inf')
    sizes = (len(node) for node in index_to_node.values())
    statistic = compute_statistic(sizes)
    heights.append(height)
    statistics.append(statistic)

    return np.array(heights), np.array(statistics)

def load_statistics(edge_list_file, index_gene_file, reverse=True):
    T = load_edge_list(edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(index_gene_file)
    return compute_statistics(T, index_to_gene, reverse)

# Run script.
def run(args):
    # Load data.
    if args.verbose:
        progress('Loading data...')

    T = load_edge_list(args.edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(args.index_gene_file)

    # Compute cluster size statistic at each height in the hierarchy.
    if args.verbose:
        progress('Computing hierarchy statistic...')
    heights, statistics = compute_statistics(T, index_to_gene)

    # Save results.
    if args.verbose:
        progress('Saving results...')

    with open(args.cluster_size_file, 'w') as f:
        f.write('\n'.join('{}\t{}'.format(height, statistic) for height, statistic in zip(heights, statistics)))

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
