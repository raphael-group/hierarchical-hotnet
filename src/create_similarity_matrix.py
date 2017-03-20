#!/usr/bin/python

# Load modules.
import math, numpy as np
import sys, argparse

from common import convert_edge_list_to_adjacency_matrix, hh_similarity_matrix
from hhio import load_edge_list, save_matrix

# Parse arguments.
def get_parser():
    description = 'Construct the Hierarchical HotNet similarity matrix with the given restart probability.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-b', '--beta', type=float, required=True, help='Restart probability')
    parser.add_argument('-n', '--name', type=str, required=False, default='PPR', help='Similarity matrix name')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output filename')
    return parser

# Run script.
def run(args):
    edge_list = load_edge_list(args.edge_list_file)
    A = convert_edge_list_to_adjacency_matrix(edge_list)
    P = hh_similarity_matrix(A, args.beta)
    save_matrix(args.output_file, P, args.name)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
