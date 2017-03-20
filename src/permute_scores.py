#!/usr/bin/python

# Load modules.
import math, numpy as np
import sys, argparse, time

from hhio import load_gene_score, load_index_gene, save_gene_score

# Parse arguments.
def get_parser():
    description = 'Permute gene scores.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--gene_score_file', type=str, required=True, help='Input filename')
    parser.add_argument('-igf', '--index_gene_file', type=str, required=False, help='Index-gene filename')
    parser.add_argument('-s', '--seed', type=int, required=False, default=10*6*time.time(), help='Random seed')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output filename')
    return parser

# Run script.
def run(args):
    # Load edges.
    gene_to_score = load_gene_score(args.gene_score_file)
    score_genes = set(gene_to_score)
    if args.index_gene_file:
        index_to_gene, gene_to_index = load_index_gene(args.index_gene_file)
        network_genes = set(gene_to_index)
    else:
        network_genes = score_genes

    # Permute scores.
    genes = sorted(score_genes)
    scores = np.array([gene_to_score[gene] for gene in genes])

    np.random.seed(args.seed)
    permute_indices = [i for i, gene in enumerate(genes) if gene in network_genes]   # Only permute scores for genes in network if given.
    scores[permute_indices] = np.random.permutation(scores[permute_indices])

    gene_to_score = dict((gene, score) for gene, score in zip(genes, scores))

    # Save permuted_scores.
    save_gene_score(args.output_file, gene_to_score)

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
