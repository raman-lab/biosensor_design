#!/usr/bin/env python

import argparse
import collections
import operator
import itertools
import math
import matplotlib.pyplot as plt
import sys

from matplotlib.backends.backend_pdf import PdfPages


def generate_sets(fastas, n):
    """generate sets of integers corresponding to seq positions. check input conditions"""
    for fasta in fastas:
        with open(fasta, 'r') as f:
            f.next()
            seq_length = len(f.next().strip())
            if not all(len(line) == seq_length + 1 for line in f if not line.startswith('>')):
                raise Exception('Error: all sequences must be same length')
    combinations = itertools.combinations(range(seq_length), n)
    return combinations


def get_subsets(fastas, index_set):
    """generates the subset of amino acids specified by the set of indices for each sequence"""
    residue_sets = []
    for fasta in fastas:
        with open(fasta, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    residue_list = []
                    for index in index_set:
                        residue_list.append(line[index])
                    residue_sets.append(''.join(residue_list))
    return residue_sets


def compute_entropy(parameters, total_count):
    """explicitly computes the entropy of frequency and entropy of sequence for a given set of positions"""
    H_k = 0
    H_s = 0
    for frequency, multiplicity in parameters.items():
        prob_km = (frequency * multiplicity) / float(total_count)
        h_k_i = prob_km * math.log(prob_km, 2)
        h_s_i = prob_km * math.log(frequency / float(total_count), 2)
        H_k += h_k_i
        H_s += h_s_i
    return -H_s, -H_k


def apply_cutoff(entropy_values, cutoff):
    """filter values such that only values above the cutoff percentile are retained"""
    number_kept = int(cutoff * len(entropy_values))
    filtered_values = collections.OrderedDict(sorted(entropy_values.items(), key=operator.itemgetter(1))[:number_kept])
    return filtered_values


def generate_plot(filtered_entropy_values):
    """scatter plot of H(K) vs H(s)"""
    almost_gray = '#808080'
    almost_black = '#262626'
    h_s, h_k = zip(*filtered_entropy_values.values())
    with PdfPages('entropy_freq_dist.pdf') as pdf:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.scatter(h_s, h_k, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                    label='before reversion')

        plt.xlim(0, max(h_s)+1)
        plt.ylim(0, max(h_k)+1)
        plt.xlabel('H[s]', fontsize=20)
        plt.ylabel('H[K]', fontsize=20)

        spines_to_remove = ['top', 'right']
        for spine in spines_to_remove:
            ax1.spines[spine].set_visible(False)
        ax1.xaxis.set_ticks_position('none')
        ax1.yaxis.set_ticks_position('none')
        spines_to_keep = ['bottom', 'left']
        for spine in spines_to_keep:
            ax1.spines[spine].set_linewidth(0.5)
            ax1.spines[spine].set_color(almost_black)
        ax1.xaxis.label.set_color(almost_black)
        ax1.yaxis.label.set_color(almost_black)

        pdf.savefig(fig)
        plt.close()


def write_output(filtered_entropy_values, fastas, n):
    """output fasta file of the form:
            >sequence name     H(s)    H(K)
            sequence subset"""
    # for index_set, entropy_tuple in filtered_entropy_values.items():
    index_set, entropy_tuple = filtered_entropy_values.items()[0]
    print index_set
    for fasta in fastas:
        output_file = '{0}_entropy_n{1}.fasta'.format(fasta.split('.')[0], n)
        with open(output_file, 'w') as o:
            with open(fasta, 'r') as f:
                for identifier, sequence in itertools.izip_longest(f, f, fillvalue=None):
                    sequence_list = list(sequence.strip())
                    for index, residue in enumerate(sequence_list):
                        if index not in index_set:
                            sequence_list[index] = '-'
                    o.write('{0}\t {1}\t {2}\n'.format(identifier.strip(), entropy_tuple[0], entropy_tuple[1]))
                    o.write('{0}\n'.format(''.join(sequence_list)))
    x = filtered_entropy_values.values()[0][0]
    y = filtered_entropy_values.values()[0][1]
    sys.stdout.write('{0}\t{1}\n'.format(x, y))


def cvs(fastas, n, cutoff, plot, greedy):
    """control flow fxn. calls child fxn. stores number of sequences M, dictionary with each unique subset as key
    and its frequency as the value, dictionary with frequency as key and multiplicity as value."""
    set_space = generate_sets(fastas, n)
    entropy_values = {}
    for index_set in set_space:
        subset_list = get_subsets(fastas, index_set)
        subset_dict = {}
        for residue_subset in subset_list:
            if residue_subset in subset_dict.keys():
                subset_dict[residue_subset] += 1
            else:
                subset_dict[residue_subset] = 1
        parameters = {}
        for frequency in subset_dict.values():
            if frequency in parameters.keys():
                parameters[frequency] += 1
            else:
                parameters[frequency] = 1
        entropy_seq, entropy_freq = compute_entropy(parameters, len(subset_list))
        entropy_values[index_set] = (entropy_seq, entropy_freq)
    filtered_entropy_values = apply_cutoff(entropy_values, cutoff)
    generate_plot(entropy_values)
    write_output(filtered_entropy_values, fastas, n)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to select maximally informative subsets of amino acid residues from protein sequences
        based on entropy of the frequency distribution and sequence of the subset."""
    )
    parser.add_argument('-g', '--greedy', action='store_true', default=False,
                        help='use greedy search algorithm with 20L iterations. '
                             'default search algorithm is combinatorial and deterministic')
    parser.add_argument('-p', '--plot', action='store_true', default=False,
                        help='generate plot of entropy functions')
    parser.add_argument('-c', '--cutoff', type=float, default=0.1,
                        help='fraction of top scoring subsets to output. should be a float from 0 to 1')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fastas', nargs='*', required=True,
                          help='fasta(s) containing protein sequences. sequences must be same length.')
    required.add_argument('-n', type=int, required=True,
                          help='size of desired subset')
    args = parser.parse_args()
    cvs(args.fastas, args.n, args.cutoff, args.plot, args.greedy)
