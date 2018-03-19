#!/usr/bin/env python
import argparse

import sys

from color_palettes import palette
import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd


def parse_fasta_to_tag_seq_dict(fasta_file):
    tag_sequence_dict = {}
    with open(fasta_file, 'r') as f:
        for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
            tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
            sequence = sequence.rstrip()
            tag_sequence_dict[tag] = sequence
    return tag_sequence_dict


def dataframe_from_file(data_file):
    columns = ['repressed', 'induced', 'free', 'fold_induction', 'fold_repression']
    with open(data_file, 'r') as f:
        dataframe = pd.read_table(f, sep='\s+', header=None, skiprows=2, index_col=0)
        dataframe.columns = columns
        dataframe.index.name = 'sequence'
    return dataframe


def histogram_from_dict(label_data_dict, num_bins, name):
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(label_data_dict)]
    with PdfPages(name) as pdf:
        fig = plt.figure()
        for k, key in enumerate(label_data_dict.iterkeys()):
            bins = np.linspace(min(label_data_dict[key]), max(label_data_dict[key]), num=num_bins)
            ax1 = fig.add_subplot(111)
            ax1.hist(label_data_dict[key], bins, alpha=0.5, color=color_set[k], label='{0}'.format(key))

        legend = ax1.legend(loc='best', scatterpoints=1, framealpha=0.5)
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

        plt.title(name, fontsize=20, y=1.02)
        plt.xlabel('Values', fontsize=20)
        plt.ylabel('Frequency', fontsize=20)
        pdf.savefig()
        plt.close()


def dna_shape_discrete(data_file, fasta_file, shape_files):
    tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta_file)
    dataframe = dataframe_from_file(data_file)
    labels = ['low', 'mid', 'high']
    dataframe['quantile'] = pd.qcut(dataframe['fold_induction'], len(labels), labels=labels)
    # sys.stdout.write(dataframe.to_string())
    for shape_file in shape_files:
        quantile_shape_dict = {}
        for label in labels:
            quantile_shape_dict[label] = []
        split_name = shape_file.split('.')
        name = '.'.join([split_name[0], split_name[-1], 'pdf'])
        with open(shape_file, 'r') as f:
            for tag, data_string in itertools.izip_longest(f, f, fillvalue=None):
                tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
                string_list = [x for x in data_string.rstrip().split(',') if x != 'NA']
                data_list = map(float, string_list)
                sequence = tag_seq_dict[tag]
                quantile = dataframe.loc[sequence]['quantile']
                quantile_shape_dict[quantile].extend(data_list)

        histogram_from_dict(quantile_shape_dict, 20, name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""script to run DNAshape parameter calculations on sequences on all 
    sequences in sequence-induction datafile from genotype_fluorescence.py. input includes sequence-induction data 
    file""")

    required = parser.add_argument_group('required arguments')
    required.add_argument('-d', '--data_file', required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    parser.add_argument('-f', '--fasta', required=True,
                        help='the fasta file containing sequences appearing in the data file')
    parser.add_argument('-s', '--shape_files', nargs='*', required=True,
                        help='shape files from DNAshape prediction. for each shape file, plots will be made')
    args = parser.parse_args()
    dna_shape_discrete(args.data_file, args.fasta, args.shape_files)
