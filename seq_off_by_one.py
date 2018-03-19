#!/usr/bin/env python
import argparse
import itertools
import pandas as pd
import sys

from induction_plotly import dataframe_from_file


def main(data_file, bins):
    allowable_lengths = [16, 17, 18, 19]
    dataframe = dataframe_from_file(data_file)
    seq_list = list(dataframe.index)
    seq_length_list = []
    for x in seq_list:
        if type(x) is str:
            seq_length_list.append(len(x))
        else:
            seq_length_list.append(0)
    dataframe['length'] = pd.Series(seq_length_list, index=dataframe.index)
    dataframe = dataframe[dataframe['length'].isin(allowable_lengths)]
    dataframe.dropna(subset=['fold_induction'], how='any', inplace=True)

    if not bins:
        labels = ['low', 'mid', 'high']
        dataframe['bins'] = pd.qcut(dataframe['fold_induction'], len(labels), labels=labels)
    else:
        dataframe['bins'] = pd.cut(dataframe['fold_induction'], bins)
    seq_list = list(dataframe.index)

    # sys.stdout.write('graph [\n')
    output_nodes = ['name\tbin\n']
    for seq in seq_list:
        bin = dataframe.loc[seq]['bins']
        output_line = '{0}\t{1}\n'.format(seq, bin)
        output_nodes.append(output_line)
    with open('nodes.txt', 'w') as o:
        o.writelines(output_nodes)
        # commented out lines are for gml formatting, which for now isnt being used
        # sys.stdout.write('\tnode [\n')
        # sys.stdout.write('\t\tid {0}\n'.format(s))
        # sys.stdout.write('\t\tlabel "{0}"\n'.format(dataframe.loc[seq]['bins']))
        # sys.stdout.write('\t]\n')

    sys.stdout.write('source\ttarget\tinteraction\tsource_attr\ttarget_attr\tdata\n')
    for seq1, seq2 in itertools.combinations(seq_list, 2):
        diff_count = sum(c1 != c2 for c1, c2 in itertools.izip_longest(seq1, seq2))
        if diff_count == 1:
            induction_list = [dataframe.loc[seq1]['fold_induction'], dataframe.loc[seq2]['fold_induction']]
            ffi = max(induction_list) / min(induction_list)
            # sys.stdout.write('\tedge [\n')
            # sys.stdout.write('\t\tlabel "{0}"\n'.format(ffi))
            # sys.stdout.write('\t\tsource {0}\n'.format(seq_list.index(seq1)))
            # sys.stdout.write('\t\ttarget {0}\n'.format(seq_list.index(seq2)))
            # sys.stdout.write('\t]\n')
            source_attr = dataframe.loc[seq1]['bins']
            target_attr = dataframe.loc[seq2]['bins']
            sys.stdout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(seq1, seq2, 'pp', source_attr, target_attr, ffi))
    # sys.stdout.write(']')

    # could do something for off by two with a different connections type (ie not pp)
    # may want to do something for edge length/ordering other than ffi. one-off edges are not organized very usefully

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to find sequences that differ by one letter and
    output a table for use in cytoscape""")
    parser.add_argument('-b', '--bins', nargs='*', type=int,
                        help='space separated integer numbers for binning fold induction values. default is to use 3 '
                             'equally populated bins')
    required = parser.add_argument_group("required")
    required.add_argument('-d', '--data', required=True, help="file with column of sequences")
    args = parser.parse_args()
    main(args.data, args.bins)
