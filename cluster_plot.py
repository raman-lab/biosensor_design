#!/usr/bin/env python
import argparse
import collections
import itertools
from color_palettes import palette
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def parse_cd_hit_cluster_file(cd_hit_cluster_file):
    cluster_seq_count_dict = collections.defaultdict()
    cluster_count = 0
    with open(cd_hit_cluster_file, 'r') as f:
        for line in f:
            if line.startswith('0'):
                cluster_count += 1
                seq_id_list = []
                cluster_size = 1
                next_line = f.next()
                while next_line[0].isdigit():
                    pre_identifier = line.rstrip().split()[2]
                    seq_id = pre_identifier.split('>')[-1].split('...')[0]
                    seq_id_list.append(seq_id)
                    cluster_size += 1
                cluster_seq_count_dict[cluster_count] = (seq_id_list, cluster_size)
    return cluster_seq_count_dict


def get_data_by_column(axes, data_file, seq_count_dict):
    possible_lengths = [16, 17, 18, 19]
    indices = {}
    indices_found = False
    # length_tuple_dict = collections.OrderedDict()
    # for length in possible_lengths:
    #     length_tuple_dict[length] = []

    length_seq_induction_dict = collections.OrderedDict()
    for length in possible_lengths:
        length_seq_induction_dict[length] = collections.defaultdict(list)

    with open(data_file, 'r') as f:
        while not indices_found:
            line = f.readline().rstrip()
            if axes[1] in line:
                split_line = line.split()
                # x_index = split_line.index(axes[0]) + 1
                y_index = split_line.index(axes[1]) + 1
                # indices[axes[0]] = x_index
                indices[axes[1]] = y_index
                indices_found = True
        for line in f:
            split_line = line.rstrip().split()
            try:
                # x = float(split_line[indices[axes[0]]])
                y = float(split_line[indices[axes[1]]])
                seq = split_line[0]
                seq_len = len(seq)
            except IndexError:
                continue
            if not math.isnan(y) and seq_len in possible_lengths:
                if seq_count_dict and seq in seq_count_dict.keys():
                    cluster_number = seq_count_dict[seq]
                    length_seq_induction_dict[seq_len][cluster_number].append(y)
    return length_seq_induction_dict


def parse_fasta_to_tag_seq_dict(fasta_file):
    tag_sequence_dict = {}
    with open(fasta_file, 'r') as f:
        for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
            tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
            sequence = sequence.rstrip()
            tag_sequence_dict[tag] = sequence
    return tag_sequence_dict


def generate_plots(data_len_tuple_dict, axes, name):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(data_len_tuple_dict) - 1]
    log_axes = ['repressed', 'induced', 'free']
    normal_axes = ['fold_induction', 'fold_repression']

    with PdfPages(name) as pdf:
        for i, length in enumerate(data_len_tuple_dict.keys()):
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            x = data_len_tuple_dict[length]
            y = (max(data_len_tuple_dict[length][x]) - min(data_len_tuple_dict[length][x])) \
                / min(data_len_tuple_dict[length][x])
            ax1.scatter(x, y, c=color_set[i-1], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                        label='{0}'.format(length))

            legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
            rect = legend.get_frame()
            rect.set_linewidth(0.25)
            texts = legend.texts
            for t in texts:
                t.set_color(almost_black)

            plt.xlabel('{0}'.format(axes[0]), fontsize=20)
            plt.ylabel(axes[1], fontsize=20)

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


def plot_score_files(data_file, cluster_file, fasta_file, y_axis):
    """create axes variable and calls functions"""
    axes = ['cluster number', y_axis]
    if cluster_file and fasta_file:
        cluster_seq_count_dict = parse_cd_hit_cluster_file(cluster_file)
        tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta_file)
        seq_count_dict = {}
        for cluster_number in cluster_seq_count_dict.iterkeys():
            seq_id_list, cluster_size = cluster_seq_count_dict[cluster_number]
            if cluster_size > 1:
                for seq_id in seq_id_list:
                    seq_count_dict[tag_seq_dict[seq_id]] = cluster_number
    else:
        seq_count_dict = False

    data_len_cluster_ind_dict = get_data_by_column(axes, data_file, seq_count_dict)
    generate_plots(data_len_cluster_ind_dict, axes, '{0}.pdf'.format(data_file.split('.')[0]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-y", "--y_axis", default='fold_induction ratio',
                        help="criterion to be plotted on y-axis (default: fold_induction)")
    parser.add_argument("-n", "--name", default='induction_plot.pdf',
                        help='name of output pdf (default: induction_plot.pdf NOW BEING IGNORED')

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    parser.add_argument('-c', '--cluster', required=True,
                        help='a sorted cluster file from cd-hit')
    parser.add_argument('-f', '--fasta', required=True,
                        help='the fasta file containing sequences appearing in the sorted cluster file')
    args = parser.parse_args()
    plot_score_files(args.data_file, args.cluster, args.fasta, args.y_axis)