#!/usr/bin/python

# Load modules.
import numpy as np
import sys, argparse
import multiprocessing as mp

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

from hhio import progress

# Parse arguments.
def get_parser():
    description = 'Plot size of clusters for observed and permuted data.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-osf', '--observed_size_file', type=str, required=True, help='Observed cluster size filename')
    parser.add_argument('-psf', '--permuted_size_files', type=str, required=True, nargs='*', help='Permuted cluster size filenames')
    parser.add_argument('-minsf', '--min_size_file', type=str, required=True, help='Minimum cluster size filename')
    parser.add_argument('-esf', '--expected_size_file', type=str, required=True, help='Expected cluster size filename')
    parser.add_argument('-maxsf', '--max_size_file', type=str, required=True, help='Maximum cluster size filename')
    parser.add_argument('-ch', '--cut_height', type=float, required=False, help='Cut height')
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1, help='Number of cores')
    parser.add_argument('-l', '--label', type=str, required=False, nargs='*', default='', help='Plot label')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose')
    return parser

# Define functions.
def load_heights_sizes(filename):
    heights = list()
    sizes = list()
    with open(filename, 'r') as f:
        for l in f:
            arrs = l.rstrip().split()
            height = float(arrs[0])
            size = float(arrs[1])
            heights.append(height)
            sizes.append(size)
    return np.array(heights), np.array(sizes)

# Run script.
def run(args):
    # Load data.
    if args.verbose:
        progress('Loading data...')

    observed_heights, observed_sizes = load_heights_sizes(args.observed_size_file)

    map_input = args.permuted_size_files

    if args.num_cores!=1:
        pool = mp.Pool(None if args.num_cores==-1 else args.num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    map_output = map_fn(load_heights_sizes, map_input)

    if args.num_cores!=1:
        pool.close()
        pool.join()

    permuted_heights_collection, permuted_sizes_collection = zip(*map_output)

    min_heights, min_sizes = load_heights_sizes(args.min_size_file)
    expected_heights, expected_sizes = load_heights_sizes(args.expected_size_file)
    max_heights, max_sizes = load_heights_sizes(args.max_size_file)

    # Plot results.
    if args.verbose:
        progress('Plotting cluster sizes...')

    ### Restrict to finite values.
    observed_heights = np.clip(observed_heights, 1e-15, float('inf'))
    min_heights = np.clip(min_heights, 1e-15, float('inf'))
    expected_heights = np.clip(expected_heights, 1e-15, float('inf'))
    max_heights = np.clip(max_heights, 1e-15, float('inf'))
    permuted_heights_collection = [np.clip(permuted_heights, 1e-15, float('inf')) for permuted_heights in permuted_heights_collection]

    ### Define colors.
    observed_color = (0.8, 0.0, 0.0)
    permuted_color = (0.0, 0.0, 0.8)
    background_permuted_color = (0.5, 0.5, 0.8)
    alpha = 0.2

    ### Plot sizes.
    plt.figure(figsize=(5, 5))

    plt.step(1.0/observed_heights, observed_sizes, where='post', c=observed_color, linewidth=2, zorder=5, label='Observed sizes')
    plt.step(1.0/expected_heights, expected_sizes, where='post', c=permuted_color, linewidth=2, zorder=4, label='Expected sizes')

    if args.cut_height:
        i = max(k for k, height in enumerate(observed_heights) if height>=args.cut_height)
        j = max(k for k, height in enumerate(expected_heights) if height>=args.cut_height)

        plt.plot((1.0/args.cut_height, 1.0/args.cut_height), (observed_sizes[i], expected_sizes[j]), c='k', linewidth=2, alpha=0.75, zorder=6, label='Chosen cut')

    plt.step(1.0/min_heights, min_sizes, where='post', c=background_permuted_color, linewidth=1, linestyle='dotted', zorder=3, label='Permuted sizes (minimum)')

    plt.step([float('nan')], [float('nan')], where='post', c=background_permuted_color, linewidth=1, zorder=1, label='Permuted sizes (all)')
    for i, (permuted_heights, permuted_sizes) in enumerate(zip(permuted_heights_collection, permuted_sizes_collection)):
        plt.step(1.0/permuted_heights, permuted_sizes, where='post', c=background_permuted_color, linewidth=0.5, alpha=alpha, zorder=1)

    plt.step(1.0/max_heights, max_sizes, where='post', c=background_permuted_color, linewidth=1, linestyle='dashed', zorder=3, label='Permuted sizes (maximum)')

    ### Set plot properties.
    plt.xlim(0.8*np.min(1.0/observed_heights[1:-1]), 1.2*np.max(1.0/observed_heights[1:-1]))
    plt.ylim(0.8, 1.2*np.max(observed_sizes))
    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'Cut height $\delta$ ($1/\delta$)')
    plt.ylabel(r'Largest cluster size at $\delta$')
    if args.label:
        plt.title(r'Cluster sizes across hierarchy cuts for' + '\n' + r'{}'.format(' '.join(args.label)))

    ax = plt.gca()
    ax.set_facecolor('white')
    plt.setp(ax.spines.values(), color='#555555')
    plt.grid(color='#555555', linestyle='dotted', alpha=0.25)

    legend = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    frame = legend.get_frame()
    frame.set_alpha(0.0)
    plt.tight_layout()

    # Save plot.
    if args.verbose:
        progress('Saving plot...')

    plt.savefig(args.output_file, bbox_inches='tight')
    plt.close()

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
