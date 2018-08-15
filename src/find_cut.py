#!/usr/bin/python

# Load modules.
import math, numpy as np
import sys, argparse

# Parse arguments.
def get_parser():
    description = 'Find cut of hierarchy.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-osf', '--observed_size_file', type=str, required=True, help='Observed height-size filename')
    parser.add_argument('-esf', '--expected_size_file', type=str, required=True, help='Expected height-size filename')
    parser.add_argument('-lbf', '--lower_bound_fraction', type=float, required=False, default=0.001, help='Lower bound fraction of nodes for cuts')
    parser.add_argument('-ubf', '--upper_bound_fraction', type=float, required=False, default=1.0, help='Upper bound fraction of nodes for cuts')
    parser.add_argument('-hf', '--height_file', type=str, required=False, help='Height file')
    parser.add_argument('-rf', '--ratio_file', type=str, required=False, help='Ratio file')
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
    expected_heights, expected_sizes = load_heights_sizes(args.expected_size_file)

    # Find cut of hierarchy.
    if args.verbose:
        progress('Finding cut of hierarchy...')

    # Find ratio at each height.
    max_observed_size = np.max(observed_sizes)
    lower_size_bound = int(math.floor(args.lower_bound_fraction*max_observed_size))
    upper_size_bound = int(math.ceil(args.upper_bound_fraction*max_observed_size))

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

    # Save results.
    if args.verbose:
        progress('Saving results...')

    if args.height_file:
        with open(args.height_file, 'w') as f:
            f.write(str(cut_height))

    if args.ratio_file:
        with open(args.ratio_file, 'w') as f:
            f.write(str(cut_ratio))

    if args.verbose:
        progress()

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
