#!/usr/bin/python

# Load modules.
import math, numpy as np, networkx as nx
import sys, argparse
from collections import defaultdict

from hhio import load_gene_score, load_index_gene, load_edge_list

# Parse arguments.
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-igf', '--index_gene_file', type=str, required=True, help='Index-gene filename')
    parser.add_argument('-elf', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-gsf', '--gene_score_file', type=str, required=True, help='Gene score filename')
    parser.add_argument('-ms', '--min_size', type=int, required=False, help='Minimum size')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output filename')
    return parser

# Run script.
def run(args):
    # Load data.
    index_to_gene, gene_to_index = load_index_gene(args.index_gene_file)
    edge_list = load_edge_list(args.edge_list_file, index_to_gene, unweighted=True)
    gene_to_score = load_gene_score(args.gene_score_file)

    # Find degrees of of network genes in subgraph induced by scored genes.
    G = nx.Graph()
    G.add_edges_from(edge_list)
    G = G.subgraph(gene_to_score)
    G = max(nx.connected_component_subgraphs(G), key=len)
    common_genes = set(G.nodes())

    degree_to_nodes = defaultdict(set)
    for node in common_genes:
        degree = G.degree(node)
        degree_to_nodes[degree].add(node)
    distinct_degrees = set(degree_to_nodes)

    # Bin genes by degree.
    if args.min_size is not None:
        min_size = args.min_size
    else:
        min_size = len(common_genes)

    bins = list()
    current_bin = list()
    for degree in sorted(distinct_degrees, reverse=True):
        current_bin += sorted(degree_to_nodes[degree])
        if len(current_bin)>=min_size:
            bins.append(current_bin)
            current_bin = list()
    if bins:
        bins[-1] += current_bin
    elif current_bin:
        bins.append(current_bin)

    # Save degree bins.
    with open(args.output_file, 'w') as f:
        f.write('\n'.join('\t'.join(current_bin) for current_bin in bins))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
