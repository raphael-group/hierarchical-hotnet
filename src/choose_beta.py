#!/usr/bin/python

# Load modules.
import numpy as np, scipy as sp, scipy.optimize
import sys, argparse

from common import convert_edge_list_to_adjacency_matrix, hh_similarity_matrix
from hhio import load_edge_list

# Parse arguments.
def get_parser():
    description = 'Choose beta to balance the distribution of random walkers between neighbors and non-neighbors.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-d', '--digits', type=int, required=False, default=2, help='Number of digits in beta')
    parser.add_argument('-t', '--threshold', type=float, required=False, default=1.0, help='Threshold for edge weights')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output filename')
    return parser

# Define functions.
def difference(A, beta, threshold):
    P = hh_similarity_matrix(A, beta)
    np.fill_diagonal(P, 0)

    r = np.sum(P[np.where(A>=threshold)])
    s = np.sum(P[np.where(A<threshold)])
    return r-s

def balanced_beta(A, threshold, digits):
    return sp.optimize.ridder(lambda beta: difference(A, beta, threshold), a=0.1**digits, b=1.0-0.1**digits, xtol=0.1**(digits+1))

# Run script.
def run(args):
    # Load graph.
    edge_list = load_edge_list(args.edge_list_file)
    A = convert_edge_list_to_adjacency_matrix(edge_list)

    # Compute beta.
    beta = balanced_beta(A, args.threshold, args.digits)

    # Save beta.
    fmt = '{:.' + str(args.digits) + 'f}'
    with open(args.output_file, 'w') as f:
        f.write(fmt.format(beta))

if __name__ == '__main__':
    run(get_parser().parse_args(sys.argv[1:]))
