#!/usr/bin/python

# Load modules.
import math, numpy as np, networkx as nx
import sys, argparse
import multiprocessing as mp

from hierarchical_clustering import find_height_to_sizes, find_cut
from common import combined_similarity_matrix
from hhio import load_index_gene, load_weighted_edge_list, progress

# Parse arguments.
def get_parser():
    description = 'Compute test statistic on hierarchical clusterings.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-oelf', '--observed_edge_list_file', type=str, required=True, help='Observed hierarchy edge list filename')
    parser.add_argument('-oigf', '--observed_index_gene_file', type=str, required=True, help='Observed hierarchy index-gene filename')
    parser.add_argument('-pelf', '--permuted_edge_list_files', type=str, required=True, nargs='*', help='Permuted hierarchy edge list filenames')
    parser.add_argument('-pigf', '--permuted_index_gene_files', type=str, required=True, nargs='*', help='Permuted hierarchy index-gene filenames')
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1, help='Number of cores')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    return parser

# Define functions.
def statistic(sizes):
    return max(sizes)

def compute_height_to_statistic(height_to_sizes):
    return dict((height, statistic(sizes)) for height, sizes in height_to_sizes.items())

def load_height_to_statistic(edge_list_file, index_gene_file):
    T = load_weighted_edge_list(edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(index_gene_file)
    height_to_sizes = find_height_to_sizes(T, index_to_gene)
    height_to_statistic = compute_height_to_statistic(height_to_sizes)
    return height_to_statistic

def load_height_to_statistic_wrapper(arrs):
    return load_height_to_statistic(*arrs)

def combine_sizes(height_to_size_collection):
    set_heights_collection = [set(x.keys()) for x in height_to_size_collection]
    set_heights = set.union(*set_heights_collection)
    heights_collection = [sorted(x, reverse=True) for x in set_heights_collection]
    heights_union = sorted(set_heights, reverse=True)

    indices = [0 for x in heights_collection]
    lengths = [len(x) for x in heights_collection]

    height_to_sizes = dict()
    for height in heights_union:
        sizes = []
        for i, heights in enumerate(heights_collection):
            while indices[i]+1<lengths[i] and heights[indices[i]+1]>=height:
                indices[i] += 1
            sizes.append(height_to_size_collection[i][heights[indices[i]]])
        height_to_sizes[height] = sizes

    return height_to_sizes

def find_max_size_difference(observed_height_to_size, permuted_height_to_size):
    height_to_size_collection = [observed_height_to_size, permuted_height_to_size]

    set_heights_collection = [set(x.keys()) for x in height_to_size_collection]
    set_heights = set.union(*set_heights_collection)
    heights_collection = [sorted(x, reverse=True) for x in set_heights_collection]
    heights_union = sorted(set_heights, reverse=True)

    indices = [0 for x in heights_collection]
    lengths = [len(x) for x in heights_collection]

    height_to_difference = dict()
    for height in heights_union:
        sizes = []
        for i, heights in enumerate(heights_collection):
            while indices[i]+1<lengths[i] and heights[indices[i]+1]>=height:
                indices[i] += 1
            sizes.append(height_to_size_collection[i][heights[indices[i]]])
        difference = float(sizes[0])/float(sizes[1])
        height_to_difference[height] = difference

    max_difference = max(height_to_difference.values())
    height = max(height for height, difference in height_to_difference.items() if difference==max_difference)

    return height

# Run script.
def run(args):
    # Load data.
    if args.verbose:
        progress('Loading data...')

    assert len(args.permuted_edge_list_files)==len(args.permuted_index_gene_files)

    if args.num_cores!=1:
        pool = mp.Pool(None if args.num_cores==-1 else args.num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    map_input = [(args.observed_edge_list_file, args.observed_index_gene_file)] + list(zip(args.permuted_edge_list_files, args.permuted_index_gene_files))
    map_output = list(map_fn(load_height_to_statistic_wrapper, map_input))

    observed_height_to_statistic = map_output[0]
    permuted_height_to_statistic_collection = map_output[1:]

    if args.num_cores!=1:
        pool.close()
        pool.join()

    # Process data.
    if args.verbose:
        progress('Processing data...')

    permuted_height_to_sizes = combine_sizes(permuted_height_to_statistic_collection)
    permuted_height_to_average_size = dict((height, np.mean(sizes)) for height, sizes in permuted_height_to_sizes.items())
    height = find_max_size_difference(observed_height_to_statistic, permuted_height_to_average_size)
    statistic = observed_height_to_statistic[height]

    T = load_weighted_edge_list(args.observed_edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(args.observed_index_gene_file)
    clusters = find_cut(T, index_to_gene, height)

    # Save data.
    if args.verbose:
        progress('Saving data...')

    permuted_height = min(x for x in permuted_height_to_sizes.keys() if x>=height)
    permuted_sizes = permuted_height_to_sizes[permuted_height]
    p_value = float(sum(1 for x in permuted_sizes if x>=statistic))/float(len(permuted_sizes))

    header_string = '#delta: {}\n#Statistic: {}\n#p-value: {}'.format(height, statistic, p_value)
    sorted_clusters = sorted(sorted(map(sorted, clusters)), key=len, reverse=True)
    output_string = '\n'.join('\t'.join(cluster) for cluster in sorted_clusters)

    with open(args.output_file, 'w') as f:
        f.write(header_string + '\n' + output_string)

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
