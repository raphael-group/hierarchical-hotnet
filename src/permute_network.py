#!/usr/bin/python

# Load modules.
import math, random, networkx as nx
import sys, argparse, datetime

from hhio import load_edge_list, save_edge_list

# Parse arguments.
def get_parser():
    description = 'Permute networks with the double edge swap algorithm.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-c', '--connected', action='store_true', help='Preserve connectivity of graph')
    parser.add_argument('-s', '--seed', type=tuple, required=False, default=datetime.datetime.now(), help='Random seed')
    parser.add_argument('-q', '--Q', type=float, required=False, default=100, help='Minimum of Q*|E| edge swaps')
    parser.add_argument('-o', '--permuted_edge_list_file', type=str, required=True, help='Permuted edge list filename')
    return parser

# Run script.
def run(args):
    # Load unpermuted network.
    edge_list = load_edge_list(args.edge_list_file)

    # Permute network.
    G = nx.Graph()
    G.add_edges_from(edge_list)

    random.seed(args.seed)
    minimum_swaps = int(math.ceil(args.Q*G.number_of_edges()))

    if not args.connected:
        G = nx.double_edge_swap(G, minimum_swaps, 2**30)
    else:
        # If G is not connected, then we perform the connected double edge swap algorithm on a
        # largest connected component of G.
        if not nx.is_connected(G):
            G = max(nx.connected_component_subgraphs(G), key=len)

        # The current connected double edge swap algorithm does not guarantee a minimum number of
        # successful edge swaps, so we enforce it.
        current_swaps = 0
        while current_swaps<minimum_swaps:
            remaining_swaps = max(minimum_swaps-current_swaps, 100)
            additional_swaps = nx.connected_double_edge_swap(G, remaining_swaps)
            current_swaps += additional_swaps

    permuted_edge_list = G.edges()

    # Save permuted_network.
    save_edge_list(args.permuted_edge_list_file, permuted_edge_list)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
