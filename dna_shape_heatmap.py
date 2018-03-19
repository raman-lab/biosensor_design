#!/usr/bin/env python
import argparse
from dna_shape_discrete import dataframe_from_file, parse_fasta_to_tag_seq_dict
import itertools
import pandas as pd
import plotly
import plotly.graph_objs as go


def dna_shape_heatmap(shape_files, data_file, fasta_file, allowable_lengths):
    tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta_file)
    dataframe = dataframe_from_file(data_file)
    length_data_dict = {}
    for length in allowable_lengths:
        length_data_dict[length] = {'fold_induction': [], 'position': [], 'data': []}

    for shape_file in shape_files:
        split_name = shape_file.split('.')
        base_name = '_'.join([split_name[0], split_name[-1]])
        with open(shape_file, 'r') as f:
            for tag, data_string in itertools.izip_longest(f, f, fillvalue=None):
                tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
                string_list = [x for x in data_string.rstrip().split(',') if x != 'NA']
                data_list = map(float, string_list)
                sequence = tag_seq_dict[tag]
                fi = dataframe.loc[sequence]['fold_induction']
                length_data_dict[len(sequence)]['fold_induction'].append(fi)
                length_data_dict[len(sequence)]['data'].append(data_list)
                length_data_dict[len(sequence)]['position'] = range(1, len(data_list) + 1)

        for length in allowable_lengths:
            data = [
                go.Heatmap(
                    # x=length_data_dict[length]['position'],
                    # y=length_data_dict[length]['fold_induction'],
                    z=length_data_dict[length]['data'],
                    colorscale='Viridis'
                )
            ]

            layout = go.Layout(
                title='{0}_{1}'.format(base_name, length),
                xaxis=dict(ticks=''),
                yaxis=dict(ticks='', nticks=10)
            )

            fig = go.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename='{0}_{1}.html'.format(base_name, length), auto_open=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to make heatmap of dna shape parameters from DNAshape
    calculation""")
    parser.add_argument('-a', '--allowable_lengths', default=[16, 17, 18, 19],
                        help='lengths of sequences to use. separate heatmap made for each length')

    required = parser.add_argument_group('required')
    required.add_argument('-s', '--shape_files', required=True, nargs='*',
                          help='files from DNAshape. separate heatmap made for each file')
    required.add_argument('-d', '--data_file', required=True,
                          help='file from genotype_fluorescence.py')
    required.add_argument('-f', '--fasta', required=True,
                          help='the fasta file containing sequences appearing in the data file')
    args = parser.parse_args()
    dna_shape_heatmap(args.shape_files, args.data_file, args.fasta, args.allowable_lengths)
