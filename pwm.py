#!/usr/bin/env python
"""Outputs a table of positional frequencies (actually counts) for each
position in the input FASTA files where there is variation."""
import argparse
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


possible = "ACDEFGHIKLMNPQRSTVWY"


def process_seq(seq, counts):
    for i, aa in enumerate(seq):
        if aa is 'X':
            continue
        else:
            pos = counts.setdefault(i + 1, {})
            pos[aa] = pos.get(aa, 0) + 1


def main(filename, counts):
    curr_seq = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == ">" and curr_seq:
                process_seq(''.join(curr_seq), counts)
                curr_seq = []
            elif line[0] != ">":
                curr_seq.append(line)

        if curr_seq:
            process_seq(''.join(curr_seq), counts)


def output(counts):
    sys.stdout.write("Pos\t" + '\t'.join(possible) + '\n')
    poses = counts.keys()
    poses.sort()
    for pos in poses:
        if len(counts[pos]) <= 1:
            continue
        sys.stdout.write(str(pos))
        N = sum(counts[pos].values())
        for aa in possible:
            p = float(counts[pos].get(aa, 0)) / N
            if p > 0:
                log_prob = round(math.log(20*p, 2), 2)
            else:
                log_prob = float('NaN')
            sys.stdout.write("\t" + str(log_prob))
        sys.stdout.write('\n')


def make_heat_map(counts, wt_dict, input_file, resfile):
    almost_black = '#262626'
    positions = np.array(list(counts.keys()))

    xi = np.arange(len(positions))
    yi = np.arange(len(possible))
    Xi, Yi = np.meshgrid(xi, yi)

    Z = np.zeros(Xi.shape)
    xbox = []
    ybox = []
    for i in xi:
        pos = i + 1
        N = sum(counts[pos].values())
        for j in yi:
            p = float(counts[pos].get(possible[j], 0)) / N
            if p > 0:
                log_prob = round(math.log(20 * p, 2), 2)
            else:
                log_prob = float('nan')
            Z[j, i] = log_prob
            try:
                if wt_dict[pos][possible[j]]:
                        xbox.append(i)
                        ybox.append(j)
            except KeyError:
                continue
    if resfile:
        resfile_positions = []
        with open(resfile, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line:
                    first_split = line.split()[0]
                    if first_split.isdigit():
                        resfile_positions.append(int(first_split))
        x_label = []
        for p, pos in enumerate(resfile_positions):
            if p + 1 in positions:
                x_label.append(pos)
    else:
        x_label = list(positions)
    y_label = list(possible)
    masked_Z = np.ma.array(Z, mask=np.isnan(Z))

    fig, ax = plt.subplots()
    heat_map = ax.pcolormesh(masked_Z, cmap='coolwarm')
    ax.set_xticks(xi + 0.5, minor=False)
    ax.set_yticks(yi + 0.5, minor=False)
    ax.set_xticklabels(x_label, minor=False)
    ax.set_yticklabels(y_label, minor=False)
    ax.xaxis.label.set_color(almost_black)
    # ax.set_xlim(0, 13)
    ax.yaxis.label.set_color(almost_black)
    ax.tick_params(direction='out', labelsize=10, top=False, right=False)

    for k in range(0, len(xbox)):
        ax.add_patch(Rectangle((xbox[k], ybox[k]), 1, 1, fill=False, edgecolor=almost_black, lw=2))
    plt.colorbar(heat_map)
    plt.savefig('{0}.png'.format(input_file.split('.')[0]), format='png', dpi=1000)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Makes position weight matrix for given fasta sequences")
    parser.add_argument("-m", "--heat_map", action="store_true",
                        help="invokes heat map visualisation. output is png image. "
                             "name is chosen from first input fasta")
    parser.add_argument("-s", "--silent", action="store_true",
                        help="do not write pwm to stdout")
    parser.add_argument("-r", "--resfile",
                        help='if resfile is provided, labels on heatmap will be positions from resfile')
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-w", "--wt_fasta", required=True,
                           help="fasta file for wild type protein seq")
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more inquiry fasta")

    args = parser.parse_args()
    counts = {}
    for filename in args.fasta:
        main(filename, counts)
    if not args.silent:
        output(counts)
    with open(args.wt_fasta, 'r') as f:
        f.next()
        wt_sequence = f.next()
    wt_dict = {}
    for index, character in enumerate(wt_sequence):
        if character is '-':
            continue
        else:
            wt_dict[index + 1] = {character: 1}
    if args.heat_map:
        make_heat_map(counts, wt_dict, args.fasta[0], args.resfile)
