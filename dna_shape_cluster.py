#!/usr/bin/env python
import argparse
import glob
import itertools
import pandas as pd
import subprocess

from dna_shape_lineplot import dna_shape_line_plotter


def parse_fasta_to_tag_seq_dict(fasta_file):
    tag_sequence_dict = {}
    with open(fasta_file, 'r') as f:
        for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
            tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
            sequence = sequence.rstrip()
            tag_sequence_dict[tag] = sequence
    return tag_sequence_dict


def parse_cdhit_to_cluster_tag_dict(cd_hit_cluster_file):
    cluster_dict = {}
    with open(cd_hit_cluster_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                cluster_number = int(line.rstrip().split()[-1])
                cluster_dict[cluster_number] = []
            else:
                pre_identifier = line.rstrip().split()[2]
                identifier_with_ellipsis = pre_identifier.split('>')[-1]
                identifier = identifier_with_ellipsis.split('...')[0]
                cluster_dict[cluster_number].append(identifier)
    return cluster_dict


def dataframe_from_file(data_file):
    columns = ['repressed', 'induced', 'free', 'fold_induction', 'fold_repression']
    with open(data_file, 'r') as f:
        dataframe = pd.read_table(f, sep='\s+', header=None, skiprows=2, index_col=0)
        dataframe.columns = columns
        dataframe.index.name = 'sequence'
    return dataframe


def dna_shape_cluster(data_file, cluster_file, fasta, cluster_number_list):
    tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta)
    cluster_tag_dict = parse_cdhit_to_cluster_tag_dict(cluster_file)
    dataframe = dataframe_from_file(data_file)
    protein_name = data_file.split('/')[-1].split('.')[0]
    for number in cluster_number_list:
        tags = cluster_tag_dict[number]
        fasta_lines = []
        fasta_name = '{0}_{1}'.format(protein_name, number)
        for tag in tags:
            seq = tag_seq_dict[tag]
            fi = round(dataframe.loc[seq]['fold_induction'], 1)
            tag_line = '>{0} {1}\n'.format(tag, fi)
            seq_line = '{0}\n'.format(seq)
            fasta_lines.extend([tag_line, seq_line])
        with open(fasta_name, 'w') as o:
            o.writelines(fasta_lines)
        dnashape_command = ['./prediction', '-i', '{0}'.format(fasta_name)]
        p = subprocess.Popen(dnashape_command)
        p.communicate()
        dnashape_files = glob.glob('{0}.*'.format(fasta_name))
        dna_shape_line_plotter(dnashape_files)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to run DNAshape parameter calcultions on sequences in
    cd-hit clustering output. input includes cd-hit sorted cluster file, fasta that was input into cd-hit,
    sequence-induction datafile, and a space separated list of cluster numbers to operate on""")

    required = parser.add_argument_group('required arguments')
    required.add_argument('-n', '--numbers', nargs='*', type=int, required=True,
                          help='space separated list of cluster numbers to use')
    required.add_argument("-d", "--data_file", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    parser.add_argument('-c', '--cluster', required=True,
                        help='a sorted cluster file from cd-hit')
    parser.add_argument('-f', '--fasta', required=True,
                        help='the fasta file containing sequences appearing in the sorted cluster file')
    args = parser.parse_args()
    dna_shape_cluster(args.data_file, args.cluster, args.fasta, args.numbers)