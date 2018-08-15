#!/usr/bin/python

# Load modules.
import numpy as np
import sys, argparse
import multiprocessing as mp

try:
    import fortran_module
    imported_fortran_module = True
except:
    imported_fortran_module = False

from hhio import progress

# Parse arguments.
def get_parser():
    description = 'Find min/expected/max cluster sizes across hierarchies.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--size_files', type=str, required=True, nargs='*', help='Cluster size filenames')
    parser.add_argument('-minsf', '--min_size_file', type=str, required=False, help='Minimum cluster size filename')
    parser.add_argument('-esf', '--expected_size_file', type=str, required=False, help='Expected cluster size filename')
    parser.add_argument('-maxsf', '--max_size_file', type=str, required=False, help='Maximum cluster size filename')
    parser.add_argument('-nc', '--num_cores', type=int, required=False, default=1, help='Number of cores')
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

    map_input = args.size_files

    if args.num_cores!=1:
        pool = mp.Pool(None if args.num_cores==-1 else args.num_cores)
        map_fn = pool.map
    else:
        map_fn = map

    map_output = map_fn(load_heights_sizes, map_input)

    if args.num_cores!=1:
        pool.close()
        pool.join()

    heights_collection, sizes_collection = zip(*map_output)

    # Summarize cluster size on hierarchies.
    if args.verbose:
        progress('Summarizing cluster sizes...')

    num_hierarchies = len(heights_collection)
    distinct_heights = np.unique(np.concatenate(heights_collection))[::-1]
    num_distinct_heights = len(distinct_heights)

    max_indices = np.zeros(num_hierarchies, dtype=np.int64)
    for i in range(num_hierarchies):
        max_indices[i] = len(heights_collection[i])

    if imported_fortran_module:
        max_index = np.max(max_indices)
        heights = np.zeros((num_hierarchies, max_index))
        sizes = np.zeros((num_hierarchies, max_index))
        for i in range(num_hierarchies):
            heights[i, 0:max_indices[i]] = heights_collection[i]
            sizes[i, 0:max_indices[i]] = sizes_collection[i]

        summary_sizes = fortran_module.summarize_sizes(distinct_heights, heights, sizes, max_indices)
        min_sizes = summary_sizes[:, 0]
        expected_sizes = summary_sizes[:, 1]
        max_sizes = summary_sizes[:, 2]

    else:
        cur_indices = np.zeros(num_hierarchies, dtype=np.int64)
        cur_sizes = np.zeros(num_hierarchies)
        min_sizes = np.zeros(num_distinct_heights)
        expected_sizes = np.zeros(num_distinct_heights)
        max_sizes = np.zeros(num_distinct_heights)

        for k in range(num_distinct_heights):
            distinct_height = distinct_heights[k]
            for i in range(num_hierarchies):
                while cur_indices[i]<max_indices[i]-1 and heights_collection[i][cur_indices[i]+1]>=distinct_height:
                    cur_indices[i] += 1
                cur_sizes[i] = sizes_collection[i][cur_indices[i]]
            min_sizes[k] = np.min(cur_sizes)
            expected_sizes[k] = np.mean(cur_sizes)
            max_sizes[k] = np.max(cur_sizes)

    # Save results.
    if args.verbose:
        progress('Saving results...')

    if args.min_size_file is not None:
        with open(args.min_size_file, 'w') as f:
            f.write('\n'.join('{}\t{}'.format(height, size) for height, size in zip(distinct_heights, min_sizes)))

    if args.expected_size_file is not None:
        with open(args.expected_size_file, 'w') as f:
            f.write('\n'.join('{}\t{}'.format(height, size) for height, size in zip(distinct_heights, expected_sizes)))

    if args.max_size_file is not None:
        with open(args.max_size_file, 'w') as f:
            f.write('\n'.join('{}\t{}'.format(height, size) for height, size in zip(distinct_heights, max_sizes)))

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
