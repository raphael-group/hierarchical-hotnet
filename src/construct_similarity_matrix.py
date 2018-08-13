#!/usr/bin/python

# Load modules.
import numpy as np, scipy as sp, scipy.optimize
import sys, argparse

from common import convert_edge_list_to_adjacency_matrix, hh_similarity_matrix
from hhio import load_edge_list, save_matrix

# Parse arguments.
def get_parser():
    description = 'Construct the Hierarchical HotNet similarity matrix and/or the restart probability beta.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-b', '--beta', type=float, required=False, help='Restart probability beta')
    parser.add_argument('-d', '--digits', type=int, required=False, default=2, help='Number of digits in beta')
    parser.add_argument('-t', '--threshold', type=float, required=False, default=1.0, help='Threshold for edge weights')
    parser.add_argument('-n', '--name', type=str, required=False, default='PPR', help='Similarity matrix name')
    parser.add_argument('-o', '--output_file', type=str, required=False, help='Similarity matrix output filename')
    parser.add_argument('-bof', '--beta_output_file', type=str, required=False, help='Beta output filename')
    return parser

# Define functions.
def difference(A, beta, threshold):
    P = hh_similarity_matrix(A, beta)
    np.fill_diagonal(P, 0)

    r = np.sum(P[np.where(A>=threshold)])
    s = np.sum(P[np.where(A<threshold)])
    return r-s

def balanced_beta(A, threshold, digits):
    import scipy as sp, scipy.optimize
    return sp.optimize.ridder(lambda beta: difference(A, beta, threshold), a=0.1**digits, b=1.0-0.1**digits, xtol=0.1**(digits+1))

# Run script.
def run(args):
    edge_list = load_edge_list(args.edge_list_file)
    A = convert_edge_list_to_adjacency_matrix(edge_list)

    if 0<args.beta<1:
        beta = args.beta
    elif args.beta is None:
        beta = balanced_beta(A, args.threshold, args.digits)
    else:
        raise ValueError('{} invalid; beta must satisfy 0 < beta < 1.'.format(args.beta))

    P = hh_similarity_matrix(A, beta)

    if args.output_file is not None:
        save_matrix(args.output_file, P, args.name)

    if args.beta_output_file is not None:
        fmt = '{:.' + str(args.digits) + 'f}'
        with open(args.beta_output_file, 'w') as f:
            f.write(fmt.format(beta))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
