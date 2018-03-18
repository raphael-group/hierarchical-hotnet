#!/usr/bin/python

# Load modules.
import math, numpy as np
import sys, argparse

# Parse arguments.
def get_parser():
    description = 'Compute p-value.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-osf', '--observed_statistic_file', type=str, required=True, help='Observed statistic filename')
    parser.add_argument('-psf', '--permuted_statistic_files', type=str, required=True, nargs='*', help='Permuted statistics filenames')
    parser.add_argument('-o', '--output_file', type=str, required=True)
    return parser

# Define functions.
def load_statistic(filename):
    with open(filename, 'r') as f:
        statistic = float(f.readlines()[0])
    return statistic

# Run script.
def run(args):
    observed_statistic = load_statistic(args.observed_statistic_file)
    permuted_statistics = [load_statistic(filename) for filename in args.permuted_statistic_files]

    num_extreme_statistics = sum(1 for permuted_statistic in permuted_statistics if permuted_statistic>=observed_statistic)
    num_total_statistics = len(permuted_statistics)
    p_value = float(num_extreme_statistics)/float(num_total_statistics)

    with open(args.output_file, 'w') as f:
        f.write(str(p_value))

if __name__=='__main__':
    run(get_parser().parse_args(sys.argv[1:]))
