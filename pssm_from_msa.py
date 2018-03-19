#!/usr/bin/env python
import argparse
import itertools
import numpy as np
import pandas as pd


def pssm_from_msa(msa_file):
    with open(msa_file, 'r') as f:
        first_tag = f.readline()
        first_aligned_seq = f.readline()
        alignment_length = len(first_aligned_seq.rstrip())

    position_freq_matrix = pd.DataFrame(data=0.25, index=['A', 'G', 'C', 'T'], columns=range(1, alignment_length + 1))
    with open(msa_file, 'r') as f:
        for tag, aligned_seq in itertools.izip_longest(f, f, fillvalue=None):
            for index, base in enumerate(aligned_seq):
                if base in position_freq_matrix.index:
                    position_freq_matrix.loc[base][index + 1] += 1
    column_sum = position_freq_matrix.sum(axis='index')
    position_prob_matrix = position_freq_matrix.divide(column_sum / 4)
    pssm = position_prob_matrix.applymap(np.log2)
    pssm.to_json('{0}_pssm.json'.format(msa_file.split('/')[-1].split('.')[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='script which makes a position-specific scoring matrix from an input multiple sequence alignment. '
                    'resulting pssm is output as a .json file with beginning with the name of the msa file. '
                    'msa should be in fasta format')
    required = parser.add_argument_group('required')
    required.add_argument('-m', '--msa', required=True, help='multiple sequence alignment file from which pssm is '
                                                             'made. msa should be fasta format')
    args = parser.parse_args()
    pssm_from_msa(args.msa)
