#!/usr/bin/python

# Load modules.
import math, numpy as np, networkx as nx
import sys, argparse
from collections import defaultdict
import multiprocessing as mp

from hierarchical_clustering import find_height_to_sizes, find_cut
from common import combined_similarity_matrix
from hhio import load_index_gene, load_weighted_edge_list, progress

# Parse arguments.
def get_parser():
    description = 'Find cut of hierarchy.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-oelf', '--observed_edge_list_file', type=str, required=True, help='Observed hierarchy edge list filename')
    parser.add_argument('-oigf', '--observed_index_gene_file', type=str, required=True, help='Observed hierarchy index-gene filename')
    parser.add_argument('-pelf', '--permuted_edge_list_files', type=str, required=True, nargs='*', help='Permuted hierarchy edge list filenames')
    parser.add_argument('-pigf', '--permuted_index_gene_files', type=str, required=True, nargs='*', help='Permuted hierarchy index-gene filenames')
    parser.add_argument('-st', '--significance_threshold', type=float, required=False, default=0.01, help='Significance threshold')
    parser.add_argument('-hf', '--height_file', type=str, required=False, help='Height file')
    parser.add_argument('-sf', '--statistic_file', type=str, required=False, help='Statistic file')
    parser.add_argument('-pf', '--p_value_file', type=str, required=False, help='p-value file')
    parser.add_argument('-rf', '--ratio_file', type=str, required=False, help='Ratio file')
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1, help='Number of cores')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    return parser

# Define functions.
def compute_statistic(sizes):
    return max(sizes)

def compute_heights(T, reverse=True):
    return sorted(set(height for source, target, height in T), reverse=reverse)

def load_heights(edge_list_file, reverse=True):
    T = load_weighted_edge_list(edge_list_file)
    return compute_heights(T, reverse=reverse)

def compute_statistics(T, index_to_gene, thresholds, reverse=True):
    # Associate nodes of dendrogram with indices, but actual node labels are unimportant for
    # computing statistic on sizes of cuts of dendrogram.
    T = sorted(T, key=lambda x: x[2], reverse=reverse)

    index_to_node = defaultdict(set)
    for index, gene in index_to_gene.items():
        index_to_node[index] = set([index])

    # For each interval, find the statistic on the hierarchy T, where each row of T has the form
    # (source, target, height).
    k = -1
    m = len(T)
    n = len(thresholds)
    statistics = list()

    for i in range(n):
        while k+1<m and (T[k+1][2]<=thresholds[i] if not reverse else T[k+1][2]>=thresholds[i]):
            k += 1
            index_to_node[T[k][1]] |= index_to_node[T[k][0]]
            del index_to_node[T[k][0]]

        # We can compute the statistic once per threshold as long as it is monotonically
        # non-decreasing.
        statistic = compute_statistic(map(len, index_to_node.values()))
        statistics.append(statistic)

    return statistics

def load_statistics(edge_list_file, index_gene_file, thresholds, reverse=True):
    T = load_weighted_edge_list(edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(index_gene_file)
    return compute_statistics(T, index_to_gene, thresholds, reverse)

def load_statistics_wrapper(arrs):
    return load_statistics(*arrs)

# Run script.
def run(args):
    assert len(args.permuted_edge_list_files)==len(args.permuted_index_gene_files)
    num_permutations = len(args.permuted_edge_list_files)

    # Compute statistic on hierarchy.
    if args.verbose:
        progress('Computing statistic on hierarchy...')

    heights = load_heights(args.observed_edge_list_file)
    observed_input = [(args.observed_edge_list_file, args.observed_index_gene_file, heights)]
    permuted_input = [(edge_list_file, index_gene_file, heights) for edge_list_file, index_gene_file in zip(args.permuted_edge_list_files, args.permuted_index_gene_files)]

    if args.num_cores!=1:
        pool = mp.Pool(None if args.num_cores==-1 else args.num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    map_input = observed_input + permuted_input
    map_output = map_fn(load_statistics_wrapper, map_input)

    if args.num_cores!=1:
        pool.close()
        pool.join()

    observed_statistics = map_output[0]
    permuted_statistics_collection = map_output[1:]

    # Find cut of hierarchy.
    if args.verbose:
        progress('Finding cut of hierarchy...')

    # Find p-values and ratio of statistics at each height.
    n = len(heights)
    height_to_statistic = dict()
    height_to_p_value = dict()
    height_to_ratio = dict()

    for i in range(n):
        height = heights[i]
        observed_statistic = observed_statistics[i]
        permuted_statistics = [statistics[i] for statistics in permuted_statistics_collection]
        p_value = float(sum(1 for permuted_statistic in permuted_statistics if permuted_statistic>=observed_statistic))/float(num_permutations)
        ratio = observed_statistic/np.mean(permuted_statistics)

        height_to_statistic[height] = observed_statistic
        height_to_p_value[height] = p_value
        height_to_ratio[height] = ratio

    # Find the largest significant height that achieves the largest ratio of observed to mean permuted statistics.
    significant_heights = set(height for height, p_value in height_to_p_value.items() if p_value<args.significance_threshold)
    maximum_ratio = max(height_to_ratio[height] for height in significant_heights)

    cut_height = max(height for height in significant_heights if height_to_ratio[height]==maximum_ratio)
    cut_statistic = height_to_statistic[cut_height]
    cut_p_value = height_to_p_value[cut_height]
    cut_ratio = height_to_ratio[cut_height]

    # Save data.
    if args.verbose:
        progress('Saving data...')

    if args.height_file:
        with open(args.height_file, 'w') as f:
            f.write(str(cut_height))

    if args.statistic_file:
        with open(args.statistic_file, 'w') as f:
            f.write(str(cut_statistic))

    if args.p_value_file:
        with open(args.p_value_file, 'w') as f:
            f.write(str(cut_p_value))

    if args.ratio_file:
        with open(args.ratio_file, 'w') as f:
            f.write(str(cut_ratio))

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
