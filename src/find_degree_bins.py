#!/usr/bin/python

# Load modules.
import math, numpy as np, networkx as nx
import sys, argparse, time
from collections import defaultdict

from hhio import load_gene_score, load_index_gene, load_edge_list

# Parse arguments.
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-igf', '--index_gene_file', type=str, required=True, help='Index-gene filename')
    parser.add_argument('-elf', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-ms', '--min_size', type=int, required=True, help='Minimum size')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output filename')
    return parser

# Run script.
def run(args):
    # Load data.
    index_to_gene, gene_to_index = load_index_gene(args.index_gene_file)
    edge_list = load_edge_list(args.edge_list_file, index_to_gene)

    # Bin genes by degree.
    G = nx.Graph()
    G.add_edges_from(edge_list)

    degree_to_nodes = defaultdict(set)
    for node in G.nodes():
        degree = G.degree(node)
        degree_to_nodes[degree].add(node)

    distinct_degrees = set(degree_to_nodes)
    degree_bins = list()
    current_degree_bin = list()
    for degree in sorted(distinct_degrees, reverse=True):
        current_degree_bin += sorted(degree_to_nodes[degree])
        if len(current_degree_bin)>args.min_size:
            degree_bins.append(current_degree_bin)
            current_degree_bin = list()
    if degree_bins:
        degree_bins[-1] += current_degree_bin
    else:
        degree_bins.append(current_degree_bin)

    # Save results.
    output_strings = ['# Rank\tMinimum degree\tMaximum degree\tNumber of genes in bin\tGenes in bin']
    for i, degree_bin in enumerate(degree_bins):
        min_degree = min(G.degree(node) for node in degree_bin)
        max_degree = max(G.degree(node) for node in degree_bin)
        output_strings.append('{}\t{}\t{}\t{}\t{}'.format(i+1, min_degree, max_degree, len(degree_bin), '|'.join(map(str, degree_bin))))

    with open(args.output_file, 'w') as f:
        f.write('\n'.join(output_strings))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
