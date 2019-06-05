#!/usr/bin/python

# Load modules.
import numpy as np, scipy as sp, scipy.optimize
import sys, argparse

from common import hh_similarity_matrix
from hhio import load_edge_list, save_matrix

# Parse arguments.
def get_parser():
    description = 'Construct the Hierarchical HotNet similarity matrix and/or the restart probability beta.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--edge_list_file', type=str, required=True, help='Edge list filename')
    parser.add_argument('-d', '--directed', action='store_true', help='Is directed')
    parser.add_argument('-b', '--beta', type=float, required=False, help='Restart probability beta')
    parser.add_argument('-nd', '--num_digits', type=int, required=False, default=2, help='Number of digits in beta')
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
    try:
        return sp.optimize.ridder(lambda beta: difference(A, beta, threshold), a=0.1**digits, b=1.0-0.1**digits, xtol=0.1**(digits+1))
    except:
        return 0.5

# Run script.
def run(args):
    # Load edge list.
    edge_list = load_edge_list(args.edge_list_file)

    # Construct adjacency matrix.
    k = min(min(edge[:2]) for edge in edge_list)
    l = max(max(edge[:2]) for edge in edge_list)

    A = np.zeros((l-k+1, l-k+1), dtype=np.float64)
    if args.directed:
        for i, j, weight in edge_list:
            A[j-k, i-k] = weight
    else:
        for i, j, weight in edge_list:
            A[i-k, j-k] = A[j-k, i-k] = weight

    # Choose beta.
    if args.beta is None:
        beta = balanced_beta(A, args.threshold, args.num_digits)
    elif 0<args.beta<1:
        beta = args.beta
    else:
        raise ValueError('{} invalid; beta must satisfy 0 < beta < 1.'.format(args.beta))

    # Construct Hierarchical HotNet similarity matrix.
    P = hh_similarity_matrix(A, beta)

    # Save results.
    if args.output_file is not None:
        save_matrix(args.output_file, P, args.name)

    if args.beta_output_file is not None:
        fmt = '{:.' + str(args.num_digits) + 'f}'
        with open(args.beta_output_file, 'w') as f:
            f.write(fmt.format(beta))

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
