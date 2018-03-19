#!/usr/bin/env python
import argparse
import itertools
import numpy as np
import sys

from dnacurve import CurvedDNA
import matplotlib.pyplot as plt
from palindrome_plot import palindrome_score


def regulondb_histogram(data_list, bins, name):
    almost_gray = '#808080'
    almost_black = '#262626'

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.hist(data_list, bins, alpha=0.5, color=almost_gray, label='n={0}'.format(len(data_list)))

    legend = ax1.legend(loc='best', framealpha=0.5)
    rect = legend.get_frame()
    rect.set_linewidth(0.0)
    texts = legend.texts
    for t in texts:
        t.set_color(almost_black)
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

    plt.title('{0} in regulon db'.format(name), fontsize=15, y=1.02)
    plt.xlabel('{0}'.format(name), fontsize=15)
    plt.ylabel('frequency', fontsize=15)
    plt.savefig('{0}.pdf'.format(name))
    plt.close()


def regulondb(regulondb_file, tfs, fasta):
    curve_angles = []
    palindrome_scores = []

    with open(regulondb_file, 'r') as f:
        for line in f:
            split_line = line.rstrip().split()
            if line.startswith('#') or len(split_line) < 12:
                continue
            tf = split_line[1]
            seq = split_line[11]
            confidence = split_line[-1]
            if confidence == 'Strong' and tf in tfs:
                binding_site = ''.join(c for c in seq if c.isupper())
                seq_curve_obj = CurvedDNA("TTGACA{0}TATAAT".format(binding_site), "NUCLEOSOME", name="Example",
                                          curve_window=10)
                curvature_angles = seq_curve_obj.curvature[0, :]
                curvature_angles = curvature_angles[curvature_angles.nonzero()]
                sum_curve_angle = curvature_angles.sum()
                curve_angles.append(sum_curve_angle)
                ps, second_half, rev_comp = palindrome_score(binding_site)
                palindrome_scores.append(ps)
    if fasta:
        with open(fasta, 'r') as f:
            for id, seq in itertools.izip_longest(f, f, fillvalue=None):
                seq = seq.rstrip()
                seq_curve_obj = CurvedDNA("TTGACA{0}TATAAT".format(seq), "NUCLEOSOME", name="Example",
                                          curve_window=10)
                curvature_angles = seq_curve_obj.curvature[0, :]
                curvature_angles = curvature_angles[curvature_angles.nonzero()]
                sum_curve_angle = curvature_angles.sum()
                curve_angles.append(sum_curve_angle)

                ps, second_half, rev_comp = palindrome_score(seq)
                palindrome_scores.append(ps)
    sys.stdout.write(','.join(map(str, palindrome_scores)))
    # curve_bins = np.linspace(min(curve_angles), max(curve_angles), num=15)
    # curve_name = 'curvature_angles'
    # regulondb_histogram(curve_angles, curve_bins, curve_name)
    # palindrome_bins = np.arange(min(palindrome_scores), max(palindrome_scores) + 1, 1)
    # palindrome_name = 'palindrome_length'
    # regulondb_histogram(palindrome_scores, palindrome_bins, palindrome_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='quick script to calculate curvature values for TF binding sites in '
                                                 'regulon database')
    parser.add_argument('-t', '--tfs',
                        default=['AcrR', 'AgaR', 'AraC', 'ArgR', 'CpxR', 'CynR', 'CytR', 'DcuR', 'GalR', 'NanR',
                                 'NarL', 'OmpR', 'OxyR', 'PhoB', 'PhoP', 'PurR', 'TrpR', 'TyrR', 'XylR', 'LacI'])
    parser.add_argument('-f', '--fasta', help='fasta file with sequences to calc curvature for')
    required = parser.add_argument_group('required')
    required.add_argument('-r', '--regulondb', help='regulon db tf binding site text file')
    args = parser.parse_args()
    regulondb(args.regulondb, args.tfs, args.fasta)
