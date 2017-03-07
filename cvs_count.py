#!/usr/bin/env python

import argparse
import itertools
import collections
import math
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def generate_sets(fastas, number):
    """generate sets of integers corresponding to seq positions. check input conditions"""
    for fasta in fastas:
        with open(fasta, 'r') as f:
            f.next()
            seq_length = len(f.next().strip())
            if not all(len(line) == seq_length + 1 for line in f if not line.startswith('>')):
                raise Exception('Error: all sequences must be same length')
    combinations = itertools.combinations(range(seq_length), number)
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


def extract_top_scores(score_dict, int_to_keep, invert, mean_minus_stdev):
    """filter values such that only values above the cutoff percentile are retained"""
    if int_to_keep <= 1:
        int_to_keep = int(round(len(score_dict) * int_to_keep))

    if invert:
        kept_scores = collections.OrderedDict(sorted(score_dict.items(), key=lambda (k, v): v[1])[:int_to_keep])
    elif mean_minus_stdev:
        mean_H_k = np.average(score_dict.values()[1])
        stdev_H_k = np.std(score_dict.values()[1])
        print mean_H_k
        print stdev_H_k
        print len(score_dict.values())
        for k, v in score_dict.items():
            # print v[1]
            if not ((mean_H_k - stdev_H_k) <= v[1] <= mean_H_k):
                del score_dict[k]
        kept_scores = collections.OrderedDict(sorted(score_dict.items(), key=lambda (k, v): v[1]))
        print len(kept_scores)
    else:
        kept_scores = collections.OrderedDict(sorted(score_dict.items(), reverse=True,
                                                     key=lambda (k, v): v[1])[:int_to_keep])
    return kept_scores


def get_counts(top_entropy_values):
    """count the number of times each residue appears in the top entropy sets"""
    counts = collections.Counter()
    for index_tuple in top_entropy_values.keys():
        tuple_count = collections.Counter(index_tuple)
        counts = counts + tuple_count
    return counts


def plot_bars(residue_counter, chart_name):
    """plot bar chart. for now this is missing because I am using it on a cluster without matplotlib"""
    almost_gray = '#808080'
    almost_black = '#262626'
    fig, ax1 = plt.subplots()
    positions = residue_counter.keys()
    counts = residue_counter.values()
    index = np.arange(len(positions))
    bar_width = 0.5

    rects1 = plt.bar(index, counts, bar_width, alpha=0.5, color=almost_gray)
    plt.ylabel('Count')
    plt.xticks(index + bar_width / 2, positions)
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

    plt.savefig("{0}.png".format(chart_name), format='png', dpi=1000)
    plt.close()


def write_counts(residue_counter):
    """print residues and their counts (if greater than one) to stdout"""
    sys.stdout.write('Residue\t Count\n')
    for residue, count in residue_counter.items():
        sys.stdout.write('{0}\t {1}\n'.format(residue, count))


def cvs_counter(fastas, largest, number, bar_chart, invert, mean_minus_stdev):
    """loop over subsets of size n. mostly uses counter objects to handle values"""
    residue_counter = collections.Counter()
    for n in number:
        set_space = generate_sets(fastas, n)
        entropy_values = {}
        for index_set in set_space:
            subset_list = get_subsets(fastas, index_set)
            subset_counter = collections.Counter(subset_list)
            freq_multiplicity = collections.Counter(subset_counter.values())
            entropy_seq, entropy_freq = compute_entropy(freq_multiplicity, len(subset_list))
            entropy_values[index_set] = (entropy_seq, entropy_freq)
        top_entropy_values = extract_top_scores(entropy_values, largest, invert, mean_minus_stdev)
        n_residue_counter = get_counts(top_entropy_values)
        residue_counter = residue_counter + n_residue_counter

    if bar_chart:
        plot_bars(residue_counter, bar_chart)
    write_counts(residue_counter)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to select maximally informative subsets of amino acid residues from protein sequences
        based on entropy of the frequency distribution and sequence of the subset.
        algorithm runs sweep from n=2 to n=<n>, selects <l> subsets with the largest entropy scores
        and outputs counts of residues to stdout, which can be optionally plotted as bar chart"""
    )
    parser.add_argument('-b', '--bar_chart', default=False,
                        help='name of bar chart, if desired, without file type. will be saved as png.')
    parser.add_argument('-n', '--number', nargs='*', type=int,  default=[5],
                        help='integer number of largest subset size to be sampled. '
                             'can accept multiple space delimited ints')
    parser.add_argument('-l', '--largest', type=float, default=0.1,
                        help='number specifying the number or percentage of of subsets with the largest entropy scores'
                             ' to keep. if number is greater than 1, that is how many subsets will be used. '
                             'if number is 1 or less, percentages will be used.')
    filtering_group = parser.add_mutually_exclusive_group()
    filtering_group.add_argument('-i', '--invert', action='store_true',
                                 help='keep the lowest scores. "largest" option still used to specify'
                                      ' how many sets to keep. only one of {invert, mean_minus_stdev} may be used.'
                                      'if neither are used, the default action is to keep '
                                      'highest scores')
    filtering_group.add_argument('-mms', '--mean_minus_stdev', action='store_true',
                                 help='keep scores that fall between the mean and one standard deviation below the'
                                      ' mean. mean - stdev <= x <= mean. only one of {invert, mean_minus_stdev} '
                                      'may be used. if neither are used, the default action is to keep highest scores.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fastas', nargs='*', required=True,
                          help='fasta(s) containing protein sequences. sequences must be the same length')
    args = parser.parse_args()
    cvs_counter(args.fastas, args.largest, args.number, args.bar_chart, args.invert, args.mean_minus_stdev)
