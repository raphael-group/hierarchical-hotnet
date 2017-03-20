#!/usr/bin/python

# Load modules.
import networkx as nx
import sys, argparse
from itertools import combinations
from collections import defaultdict

from hhio import load_index_gene, load_edge_list, progress

# Parse arguments.
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-cf', '--component_files', type=str, required=True, nargs='*')
    parser.add_argument('-igf', '--index_gene_files', type=str, required=True, nargs='*')
    parser.add_argument('-elf', '--edge_list_files', type=str, required=True, nargs='*')
    parser.add_argument('-n', '--networks', type=str, required=True, nargs='*')
    parser.add_argument('-s', '--scores', type=str, required=True, nargs='*')
    parser.add_argument('-t', '--threshold', type=int, required=True)
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbose', action='store_true')
    return parser

# Define functions.
def load_components(component_file):
    components = []
    with open(component_file, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.rstrip('\n').split('\t')
                component = sorted(arrs)
                components.append(component)
                break
    return components

# Run script.
def run(args):
    # Load data.
    if args.verbose:
        progress('Loading data...')

    assert len(args.component_files)==len(args.index_gene_files)==len(args.edge_list_files)==len(args.networks)==len(args.scores)

    index_to_gene_collection = dict()
    gene_to_index_collection = dict()
    edge_list_collection = dict()
    components_collection = dict()

    for network_label, score_label, index_gene_file, edge_list_file, component_file in zip(args.networks, args.scores, args.index_gene_files, args.edge_list_files, args.component_files):
        index_to_gene, gene_to_index = load_index_gene(index_gene_file)
        edge_list = set(frozenset((index_to_gene[i], index_to_gene[j])) for i, j in load_edge_list(edge_list_file))
        components = load_components(component_file)

        index_to_gene_collection[(network_label, score_label)] = index_to_gene
        gene_to_index_collection[(network_label, score_label)] = gene_to_index
        edge_list_collection[(network_label, score_label)] = edge_list
        components_collection[(network_label, score_label)] = components

    # Process data.
    if args.verbose:
        progress('Processing data...')

    edge_to_networks = defaultdict(set)
    edge_to_scores = defaultdict(set)
    edge_to_pairs = defaultdict(set)
    edge_to_tally = defaultdict(int)

    for network_label, score_label in zip(args.networks, args.scores):
        edge_list = edge_list_collection[(network_label, score_label)]
        components = components_collection[(network_label, score_label)]
        for component in components:
            for u, v in combinations(component, 2):
                edge = frozenset((u, v))
                if edge in edge_list:
                    edge_to_tally[edge] += 1

    thresholded_edges = set(edge for edge, tally in edge_to_tally.items() if tally>=args.threshold)

    G = nx.Graph()
    G.add_edges_from(thresholded_edges)
    consensus_results = sorted(sorted([sorted(x) for x in nx.connected_components(G)]), key=len, reverse=True)

    # Save data.
    if args.verbose:
        progress('Saving data...')

    output_string = '\n'.join('\t'.join(x) for x in consensus_results)
    with open(args.output_file, 'w') as f:
        f.write(output_string)

    if args.verbose:
        progress()

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
