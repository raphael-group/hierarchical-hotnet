#!/usr/bin/python

# Load modules.
import math, numpy as np, networkx as nx
import sys, argparse
from collections import defaultdict
import multiprocessing as mp

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

from hierarchical_clustering import find_height_to_sizes
from common import combined_similarity_matrix
from hhio import load_index_gene, load_weighted_edge_list, progress

# Parse arguments.
def get_parser():
    description = 'Plot test statistic on hierarchical clusterings.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-oelf', '--observed_edge_list_file', type=str, required=True, help='Observed hierarchy edge list filename')
    parser.add_argument('-oigf', '--observed_index_gene_file', type=str, required=True, help='Observed hierarchy index-gene filename')
    parser.add_argument('-pelf', '--permuted_edge_list_files', type=str, required=True, nargs='*', help='Permuted hierarchy edge list filenames')
    parser.add_argument('-pigf', '--permuted_index_gene_files', type=str, required=True, nargs='*', help='Permuted hierarchy index-gene filenames')
    parser.add_argument('-ch', '--cut_height', type=float, required=False, help='Cut height')
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1, help='Number of cores')
    parser.add_argument('-l', '--label', type=str, required=False, nargs='*', default='', help='Plot label')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file')
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

    observed_heights = np.asarray(load_heights(args.observed_edge_list_file))
    permuted_heights = observed_heights

    observed_input = [(args.observed_edge_list_file, args.observed_index_gene_file, observed_heights)]
    permuted_input = [(edge_list_file, index_gene_file, permuted_heights) for edge_list_file, index_gene_file in zip(args.permuted_edge_list_files, args.permuted_index_gene_files)]

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
    n = len(permuted_heights)
    permuted_statistics_min = list()
    permuted_statistics_mean = list()
    permuted_statistics_max = list()

    for i in range(n):
        permuted_statistics = np.array([statistics[i] for statistics in permuted_statistics_collection])
        permuted_statistics_min.append(np.min(permuted_statistics))
        permuted_statistics_mean.append(np.mean(permuted_statistics))
        permuted_statistics_max.append(np.max(permuted_statistics))

    # Plot results.
    if args.verbose:
        progress('Plotting data...')

    ### Plot sizes.
    inverse_observed_heights = 1.0/observed_heights
    inverse_permuted_heights = 1.0/permuted_heights
    num_permutations = len(permuted_statistics_collection)

    observed_color = (0.8, 0.0, 0.0)
    permuted_color = (0.0, 0.0, 0.8)
    background_permuted_color = (0.5, 0.5, 0.8)
    alpha = max(1.0/float(num_permutations), 0.1)

    plt.figure(figsize=(5, 5))

    plt.step(inverse_observed_heights, observed_statistics, c=observed_color, linewidth=2, zorder=5, label='Observed')
    plt.step(inverse_permuted_heights, permuted_statistics_max, c=permuted_color, linewidth=1, linestyle='dashed', zorder=3, label='Permuted (max.)')
    plt.step(inverse_permuted_heights, permuted_statistics_mean, c=permuted_color, linewidth=2, zorder=4, label='Permuted (mean)')
    plt.step(inverse_permuted_heights, permuted_statistics_min, c=permuted_color, linewidth=1, linestyle='dotted', zorder=3, label='Permuted (min.)')

    for statistics in permuted_statistics_collection:
        plt.step(inverse_permuted_heights, statistics, c=background_permuted_color, linewidth=0.5, alpha=alpha, zorder=1)

    ### Plot cut.
    if args.cut_height:
        observed_larger_height = min(height for height in observed_heights if height>=args.cut_height)
        permuted_larger_height = min(height for height in permuted_heights if height>=args.cut_height)
        observed_statistic = observed_statistics[list(observed_heights).index(observed_larger_height)]
        permuted_statistic_mean = permuted_statistics_mean[list(permuted_heights).index(permuted_larger_height)]

        plt.plot((1.0/args.cut_height, 1.0/args.cut_height),
            (observed_statistic, permuted_statistic_mean),
            c='k', ls='--', alpha=0.75, zorder=2)

    ### Set plot properties, save plot, etc.
    plt.xlim(0.8*np.min(inverse_observed_heights), 1.2*np.max(inverse_observed_heights))
    plt.ylim(0.8, 1.2*np.max(observed_statistics))
    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'$1/\delta$')
    plt.ylabel(r'Cluster size')
    if args.label:
        plt.title(r'Cluster size across $\delta$ thresholds for' + '\n' + r'{}'.format(' '.join(args.label)))

    ax = plt.gca()
    ax.set_facecolor('white')
    plt.setp(ax.spines.values(), color='#555555')
    plt.grid(color='#555555', linestyle='dotted', alpha=0.25)

    legend = plt.legend(loc='lower right', ncol=1)
    frame = legend.get_frame()
    frame.set_alpha(0.0)
    plt.tight_layout()

    if args.verbose:
        progress('Saving data...')

    if args.output_file.endswith('.png'):
        plt.savefig(args.output_file, dpi=600, bbox_inches='tight')
    else:
        plt.savefig(args.output_file, bbox_inches='tight')
    plt.close()

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
