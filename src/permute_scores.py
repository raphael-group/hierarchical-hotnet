#!/usr/bin/python

# Load modules.
import numpy as np
import sys, argparse

from hhio import load_gene_score, save_gene_score

# Parse arguments.
def get_parser():
    description = 'Permute gene scores.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--gene_score_file', type=str, required=True, help='Input filename')
    parser.add_argument('-bf', '--bin_file', type=str, required=False, help='Bin filename')
    parser.add_argument('-s', '--seed', type=int, required=False, help='Random seed')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output filename')
    return parser

def load_bins(filename):
    bins = list()
    with open(filename, 'r') as f:
        for l in f:
            if not l.startswith('#'):
                arrs = l.rstrip('\n').split('\t')
                current_bin = set(arrs)
                bins.append(current_bin)
    return bins

# Run script.
def run(args):
    # Load unpermuted scores.
    gene_to_score = load_gene_score(args.gene_score_file)

    if args.bin_file:
        bins = load_bins(args.bin_file)
    else:
        bins = [sorted(gene_to_score)]

    # Permute scores.
    genes = sorted(gene_to_score)
    scores = np.array([gene_to_score[gene] for gene in genes])

    if args.seed is not None:
        np.random.seed(args.seed)

    for permute_genes in bins:
        permute_indices = [index for index, gene in enumerate(genes) if gene in permute_genes]
        scores[permute_indices] = np.random.permutation(scores[permute_indices])

    gene_to_score = dict((gene, score) for gene, score in zip(genes, scores))

    # Save permuted_scores.
    save_gene_score(args.output_file, gene_to_score)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
