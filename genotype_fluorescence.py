#!/usr/bin/env python
import argparse
import collections
import numpy as np
import pandas as pd
import sys


def parse_barcode_file(barcode_file):
    barcode_dict = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            barcode, median_fluor = line.rstrip().split()
            barcode = barcode.upper()
            barcode_dict[barcode] = int(median_fluor)
    return barcode_dict


def check_barcode(full_sequence, barcode_dict):
    for barcode in barcode_dict.keys():
        if full_sequence.startswith(barcode) or full_sequence.endswith(barcode):
            return barcode
    return None


def calculate_genotype_fluorescence(sequence_files, barcode_file):
    barcode_dict = parse_barcode_file(barcode_file)
    sequence_dict = collections.defaultdict(dict)

    # variables for my own use - will be removed in future
    total_sequence_count = 0
    sequences_without_barcode_count = 0

    for seq_file in sequence_files:
        # sys.stdout.write('processing {}\n'.format(seq_file))
        with open(seq_file, 'r') as f:
            for line in f:
                full_sequence = line.rstrip()
                total_sequence_count += 1
                barcode = check_barcode(full_sequence, barcode_dict)
                if barcode:
                    sequence = full_sequence.replace(barcode, '')
                    if barcode not in sequence_dict[sequence].keys():
                        sequence_dict[sequence] = collections.Counter()
                    sequence_dict[sequence][barcode] += 1
                else:
                    sequences_without_barcode_count += 1

    count_matrix = pd.DataFrame.from_dict(sequence_dict, orient='index') # c_ij
    barcode_counts = count_matrix.sum(axis='index')
    barcode_fractions = barcode_counts / sum(barcode_counts) # f_j
    fraction_matrix = barcode_fractions * count_matrix.div(count_matrix.sum(axis='index'), axis='columns') # b_ij
    normalized_fraction_matrix = fraction_matrix.div(fraction_matrix.sum(axis='columns'), axis='index') # a_ij

    barcode_fluor_dict = dict((k, v) for k, v in barcode_dict.items() if k in normalized_fraction_matrix.columns.values)
    median_fluor = pd.Series.from_array(barcode_fluor_dict.values(), index=barcode_fluor_dict.keys()) # m_j

    log_fluor_matrix = normalized_fraction_matrix * np.log10(median_fluor)
    seq_vec = log_fluor_matrix.sum(axis='columns')
    mean_fluor = np.power(seq_vec, 10)

    sys.stdout.write('{}\n'.format(mean_fluor))
    sys.stdout.write('total seqs: {}\n'.format(total_sequence_count))
    sys.stdout.write('missed: {}\n'.format(sequences_without_barcode_count))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to calculate the mean fluorescence value for a specific genotype present in sequence data.
        input is text file(s) of sequences and a tab delimited text file of barcodes and median fluorescence values
        for each barcode."""
    )
    required = parser.add_argument_group('reguired arguments')
    required.add_argument('-s', '--sequences', nargs='*', help='text file(s) of sequences')
    required.add_argument('-b', '--barcodes', help='tab delimited text file with two columns. one column of barcodes '
                                                   'then one column of median fluorescence values for each barcode')
    args = parser.parse_args()
    calculate_genotype_fluorescence(args.sequences, args.barcodes)
