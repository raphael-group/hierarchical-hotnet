#!/usr/bin/python

# Load modules.
import numpy as np
import sys, argparse
from collections import defaultdict
import multiprocessing as mp

try:
    import fortran_module
    imported_fortran_module = True
except:
    imported_fortran_module = False

from hhio import load_index_gene, load_edge_list, progress

# Parse arguments.
def get_parser():
    description = 'Process observed and permuted hierarchies.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-oelf', '--observed_edge_list_file', type=str, required=True, help='Observed edge list filename')
    parser.add_argument('-oigf', '--observed_index_gene_file', type=str, required=True, help='Observed index-gene filename')
    parser.add_argument('-pelf', '--permuted_edge_list_files', type=str, required=True, nargs='*', help='Permuted edge list filenames')
    parser.add_argument('-pigf', '--permuted_index_gene_files', type=str, required=True, nargs='*', help='Permuted index-gene filenames')
    parser.add_argument('-lsb', '--lower_size_bound', type=float, required=False, default=10.0, help='Lower bound for cut size')
    parser.add_argument('-usb', '--upper_size_bound', type=float, required=False, default=float('inf'), help='Upper bound for cut size')
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1, help='Number of cores')
    parser.add_argument('-cf', '--cluster_file', type=str, required=False, help='Cluster filename')
    parser.add_argument('-pf', '--plot_file', type=str, required=False, help='Plot filename')
    parser.add_argument('-pl', '--plot_label', type=str, required=False, nargs='*', help='Plot label')
    parser.add_argument('-osf', '--observed_size_file', type=str, required=False, help='Observed cluster size filename')
    parser.add_argument('-esf', '--expected_size_file', type=str, required=False, help='Expected cluster size filename')
    parser.add_argument('-minsf', '--min_size_file', type=str, required=False, help='Minimum cluster size filename')
    parser.add_argument('-maxsf', '--max_size_file', type=str, required=False, help='Maximum cluster size filename')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    return parser

# Define functions.
def compute_statistic(sizes):
    return max(sizes)

def compute_statistics(T, index_to_gene, reverse=True):
    heights = list()
    statistics = list()

    # Associate nodes of dendrogram with indices.
    index_to_node = defaultdict(set)
    for index, gene in index_to_gene.items():
        index_to_node[index] = set([gene])

    # Add cluster statistic of leaf nodes.
    height = float('inf') if reverse else 0.0
    sizes = (len(node) for node in index_to_node.values())
    statistic = compute_statistic(sizes)
    heights.append(height)
    statistics.append(statistic)

    # Add cluster statistics for inner nodes.
    T = sorted(T, key=lambda x: x[2], reverse=reverse)
    m = len(T)

    for i, edge in enumerate(T):
        source, target, height = edge

        a = index_to_node[source]
        b = index_to_node[target]
        c = set(a) | set(b)

        del index_to_node[source]
        index_to_node[target] = c

        if i==m-1 or height!=T[i+1][2]:
            sizes = (len(node) for node in index_to_node.values())
            statistic = compute_statistic(sizes)
            heights.append(height)
            statistics.append(statistic)

    # Add cluster statistic for root node.
    height = 0.0 if reverse else float('inf')
    sizes = (len(node) for node in index_to_node.values())
    statistic = compute_statistic(sizes)
    heights.append(height)
    statistics.append(statistic)

    return np.array(heights), np.array(statistics)

def find_statistics(edge_list_file, index_gene_file, reverse=True):
    T = load_edge_list(edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(index_gene_file)
    return compute_statistics(T, index_to_gene, reverse)

def find_statistics_wrapper(arrs):
    return find_statistics(*arrs)

def find_cluster_sizes(observed_edge_list_file, observed_index_gene_file, permuted_edge_list_files, permuted_index_gene_files, num_cores=1):
    map_input = [(observed_edge_list_file, observed_index_gene_file)] + \
                [(permuted_edge_list_file, permuted_index_gene_file) for permuted_edge_list_file, permuted_index_gene_file in zip(permuted_edge_list_files, permuted_index_gene_files)]

    if num_cores!=1:
        pool = mp.Pool(None if num_cores==-1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    map_output = list(map_fn(find_statistics_wrapper, map_input))

    if num_cores!=1:
        pool.close()
        pool.join()

    observed_heights, observed_sizes = map_output[0]
    permuted_heights_collection, permuted_sizes_collection = zip(*map_output[1:])

    return observed_heights, observed_sizes, permuted_heights_collection, permuted_sizes_collection

def summarize_cluster_sizes(permuted_heights_collection, permuted_sizes_collection):
    # Find expected, minimum, and maximum sizes on permuted hierarchies across
    # all distinct cuts of permuted hierarchies.
    num_permutations = len(permuted_heights_collection)
    distinct_heights = np.unique(np.concatenate(permuted_heights_collection))[::-1]
    num_distinct_heights = len(distinct_heights)

    max_indices = np.zeros(num_permutations, dtype=np.int64)
    for i in range(num_permutations):
        max_indices[i] = len(permuted_heights_collection[i])

    # Use Fortran module if available for better performance.
    if imported_fortran_module:
        max_index = np.max(max_indices)
        heights = np.zeros((num_permutations, max_index))
        sizes = np.zeros((num_permutations, max_index))
        for i in range(num_permutations):
            heights[i, 0:max_indices[i]] = permuted_heights_collection[i]
            sizes[i, 0:max_indices[i]] = permuted_sizes_collection[i]

        summary_sizes = fortran_module.summarize_sizes(distinct_heights, heights, sizes, max_indices)
        min_sizes = summary_sizes[:, 0]
        expected_sizes = summary_sizes[:, 1]
        max_sizes = summary_sizes[:, 2]

    else:
        cur_indices = np.zeros(num_permutations, dtype=np.int64)
        cur_sizes = np.zeros(num_permutations)
        min_sizes = np.zeros(num_distinct_heights)
        expected_sizes = np.zeros(num_distinct_heights)
        max_sizes = np.zeros(num_distinct_heights)

        for k in range(num_distinct_heights):
            distinct_height = distinct_heights[k]
            for i in range(num_permutations):
                while cur_indices[i]<max_indices[i]-1 and permuted_heights_collection[i][cur_indices[i]+1]>=distinct_height:
                    cur_indices[i] += 1
                cur_sizes[i] = permuted_sizes_collection[i][cur_indices[i]]
            min_sizes[k] = np.min(cur_sizes)
            expected_sizes[k] = np.mean(cur_sizes)
            max_sizes[k] = np.max(cur_sizes)

    return distinct_heights, min_sizes, expected_sizes, max_sizes

def find_cut(observed_heights, observed_sizes, expected_heights, expected_sizes, lower_size_bound, upper_size_bound):
    num_observed_heights = len(observed_heights)
    num_expected_heights = len(expected_heights)
    ratios = np.zeros(num_observed_heights)

    j = 0
    for i in range(num_observed_heights):
        if lower_size_bound<=observed_sizes[i]<=upper_size_bound:
            while j<num_expected_heights-1 and expected_heights[j+1]>=observed_heights[i]:
                j += 1
            ratios[i] = float(observed_sizes[i])/float(expected_sizes[j])

    max_index = np.argmax(ratios)
    cut_height = observed_heights[max_index]
    cut_ratio = ratios[max_index]

    return cut_height, cut_ratio

def find_cut_wrapper(arrs):
    return find_cut(*arrs)

def evaluate_cut(observed_heights, observed_sizes, permuted_heights_collection, permuted_sizes_collection, expected_heights, expected_sizes, lower_size_bound, upper_size_bound, num_cores=1):
    map_input = [(observed_heights, observed_sizes, expected_heights, expected_sizes, lower_size_bound, upper_size_bound)] + \
                [(permuted_heights, permuted_sizes, expected_heights, expected_sizes, lower_size_bound, upper_size_bound) for permuted_heights, permuted_sizes in zip(permuted_heights_collection, permuted_sizes_collection)]

    if num_cores!=1:
        pool = mp.Pool(None if num_cores==-1 else num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    map_output = list(map_fn(find_cut_wrapper, map_input))

    if num_cores!=1:
        pool.close()
        pool.join()

    observed_cut_height, observed_cut_ratio = map_output[0]
    permuted_cut_height_collection, permuted_cut_ratio_collection = zip(*map_output[1:])

    expected_cut_ratio = np.mean(permuted_cut_ratio_collection)
    num_extreme_permuted_cut_ratios = sum(1 for permuted_cut_ratio in permuted_cut_ratio_collection if permuted_cut_ratio>=observed_cut_ratio)
    num_permuted_cut_ratios = len(permuted_cut_ratio_collection)

    p_value = float(num_extreme_permuted_cut_ratios)/float(num_permuted_cut_ratios)

    return observed_cut_height, observed_cut_ratio, expected_cut_ratio, p_value

def find_clusters(preordered_T, index_to_gene, threshold, reverse=True):
    # Associate nodes of dendrogram with indices.
    index_to_node = defaultdict(set)
    for index, gene in index_to_gene.items():
        index_to_node[index] = frozenset([gene])

    # Initialize clusters with leaf nodes.
    clusters = set(index_to_node.values())

    # Update clusters while ascending dendrogram.
    T = sorted(preordered_T, key=lambda x: x[2], reverse=reverse)
    m = len(T)

    for k, edge in enumerate(T):
        source, target, height = edge
        if (not reverse and height>threshold) or (reverse and height<threshold):
            break

        a = index_to_node[source]
        b = index_to_node[target]
        c = frozenset(set(a) | set(b))

        clusters.discard(a)
        clusters.discard(b)
        clusters.add(c)

        del index_to_node[source]
        index_to_node[target] = c

    return clusters

def cut_hierarchy(edge_list_file, index_gene_file, cut_height):
    T = load_edge_list(edge_list_file)
    index_to_gene, gene_to_index = load_index_gene(index_gene_file)
    clusters = find_clusters(T, index_to_gene, cut_height)

    return clusters

def plot_cluster_sizes(observed_heights, observed_sizes, min_heights, min_sizes, expected_heights, expected_sizes, max_heights, max_sizes, permuted_heights_collection, permuted_sizes_collection, cut_height, plot_label, plot_file):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')

    # Restrict to finite values.
    observed_heights = np.clip(observed_heights, 1e-15, float('inf'))
    min_heights = np.clip(min_heights, 1e-15, float('inf'))
    expected_heights = np.clip(expected_heights, 1e-15, float('inf'))
    max_heights = np.clip(max_heights, 1e-15, float('inf'))
    permuted_heights_collection = [np.clip(permuted_heights, 1e-15, float('inf')) for permuted_heights in permuted_heights_collection]

    # Define colors.
    observed_color = (0.8, 0.0, 0.0)
    permuted_color = (0.0, 0.0, 0.8)
    background_permuted_color = (0.5, 0.5, 0.8)
    alpha = 0.2

    # Plot sizes.
    plt.figure(figsize=(5, 5))

    plt.step(1.0/observed_heights, observed_sizes, where='post', c=observed_color, linewidth=2, zorder=5, label='Observed sizes')
    plt.step(1.0/expected_heights, expected_sizes, where='post', c=permuted_color, linewidth=2, zorder=4, label='Expected sizes')

    if cut_height:
        i = max(k for k, height in enumerate(observed_heights) if height>=cut_height)
        j = max(k for k, height in enumerate(expected_heights) if height>=cut_height)

        plt.plot((1.0/cut_height, 1.0/cut_height), (observed_sizes[i], expected_sizes[j]), c='k', linewidth=2, alpha=0.75, zorder=6, label='Chosen cut')

    plt.step(1.0/min_heights, min_sizes, where='post', c=background_permuted_color, linewidth=1, linestyle='dotted', zorder=3, label='Permuted sizes (minimum)')

    plt.step([float('nan')], [float('nan')], where='post', c=background_permuted_color, linewidth=1, zorder=1, label='Permuted sizes (all)')
    for i, (permuted_heights, permuted_sizes) in enumerate(zip(permuted_heights_collection, permuted_sizes_collection)):
        plt.step(1.0/permuted_heights, permuted_sizes, where='post', c=background_permuted_color, linewidth=0.5, alpha=alpha, zorder=1)

    plt.step(1.0/max_heights, max_sizes, where='post', c=background_permuted_color, linewidth=1, linestyle='dashed', zorder=3, label='Permuted sizes (maximum)')

    # Set plot properties.
    plt.xlim(0.8*np.min(1.0/observed_heights[1:-1]), 1.2*np.max(1.0/observed_heights[1:-1]))
    plt.ylim(0.8, 1.2*np.max(observed_sizes))
    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'Cut height $\delta$ ($1/\delta$)')
    plt.ylabel(r'Largest cluster size at $\delta$')
    if plot_label:
        plt.title(r'Cluster sizes across hierarchy cuts for' + '\n' + r'{}'.format(plot_label))

    ax = plt.gca()
    ax.set_facecolor('white')
    plt.setp(ax.spines.values(), color='#555555')
    plt.grid(color='#555555', linestyle='dotted', alpha=0.25)

    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    frame = legend.get_frame()
    frame.set_alpha(0.0)
    plt.tight_layout()

    # Save plot.
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close()

# Run script.
def run(args):
    # Load data.
    if args.verbose:
        progress('Loading hierarchy data...')

    observed_heights, observed_sizes, permuted_heights_collection, permuted_sizes_collection = find_cluster_sizes(args.observed_edge_list_file, args.observed_index_gene_file, args.permuted_edge_list_files, args.permuted_index_gene_files, args.num_cores)

    # Summarize hierarchy cluster sizes.
    if args.verbose:
        progress('Summarizing hierarchy cluster sizes...')

    distinct_heights, min_sizes, expected_sizes, max_sizes = summarize_cluster_sizes(permuted_heights_collection, permuted_sizes_collection)

    if args.observed_size_file is not None:
        with open(args.observed_size_file, 'w') as f:
            f.write('\n'.join('{}\t{}'.format(height, size) for height, size in zip(observed_heights, observed_sizes)))

    if args.expected_size_file is not None:
        with open(args.expected_size_file, 'w') as f:
            f.write('\n'.join('{}\t{}'.format(height, size) for height, size in zip(distinct_heights, expected_sizes)))

    if args.min_size_file is not None:
        with open(args.min_size_file, 'w') as f:
            f.write('\n'.join('{}\t{}'.format(height, size) for height, size in zip(distinct_heights, min_sizes)))

    if args.max_size_file is not None:
        with open(args.max_size_file, 'w') as f:
            f.write('\n'.join('{}\t{}'.format(height, size) for height, size in zip(distinct_heights, max_sizes)))

    # Find hierarchy cut.
    if args.verbose:
        progress('Finding hierarchy cuts...')

    observed_cut_height, observed_cut_ratio, expected_cut_ratio, p_value = evaluate_cut(observed_heights, observed_sizes, permuted_heights_collection, permuted_sizes_collection, distinct_heights, expected_sizes, args.lower_size_bound, args.upper_size_bound, args.num_cores)

    # Cut hierarchy.
    if args.verbose:
        progress('Cutting hierarchies...')

    observed_clusters = cut_hierarchy(args.observed_edge_list_file, args.observed_index_gene_file, observed_cut_height)
    observed_cut_size = compute_statistic(map(len, observed_clusters))
    expected_cut_size = observed_cut_size/observed_cut_ratio if observed_cut_ratio else float('nan')

    if args.cluster_file is not None:
        sorted_observed_clusters = sorted(sorted(map(sorted, observed_clusters)), key=len, reverse=True)
        observed_cluster_string = '\n'.join('\t'.join(cluster) for cluster in sorted_observed_clusters)
        cluster_string = '# Observed cut height: {}\n# Observed size of largest cluster at observed cut height: {}\n# Expected size of largest cluster at observed cut height: {}\n# Observed maximum ratio statistic: {:.3f}\n# Expected maximum ratio statistic: {:.3f}\n# p-value: {}\n# Clusters:\n{}'.format(observed_cut_height, observed_cut_size, expected_cut_size, observed_cut_ratio, expected_cut_ratio, p_value, observed_cluster_string)

        with open(args.cluster_file, 'w') as f:
            f.write(cluster_string)

    # Plot cluster sizes.
    if args.plot_file is not None:
        if args.verbose:
            progress('Plotting cluster sizes...')

        plot_cluster_sizes(observed_heights, observed_sizes, distinct_heights, min_sizes, distinct_heights, expected_sizes, distinct_heights, max_sizes, permuted_heights_collection, permuted_sizes_collection, observed_cut_height, ' '.join(args.plot_label), args.plot_file)

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
