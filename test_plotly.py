#!/usr/bin/env python
import argparse
import collections
from color_palettes import palette
import pandas as pd
import plotly
import math
import matplotlib
import plotly.graph_objs as go
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def generate_plots(dataframe, axes, color_scheme, name, color_dict):
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(dataframe.columns.values)]
    log_axes = ['repressed', 'induced', 'free']
    normal_axes = ['fold_induction', 'fold_repression']

    dataframe.dropna(subset=[axes], how='any', inplace=True)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    x = dataframe[axes[0]]
    y = dataframe[axes[1]]
    color_list = [color_dict[key] for key in dataframe[color_scheme]]
    ax1.scatter(x, y, c=color_list,
                marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15, label=dataframe[color_scheme])
    legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
    rect = legend.get_frame()
    rect.set_linewidth(0.25)
    texts = legend.texts
    for t in texts:
        t.set_color(almost_black)
    if axes[0] in log_axes:
        plt.xscale('log')
        ax1.set_xbound(lower=100.0)
    else:
        ax1.set_xbound(lower=0.0)
    if axes[1] in log_axes:
        plt.yscale('log')
        ax1.set_ybound(lower=100.0)
    else:
        ax1.set_ybound(lower=0.0)
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
    plot_url = py.plot_mpl(fig, filename="mpl-scatter")
    plt.close()


def get_data_by_column(axes, data_file, seq_count_dict):
    possible_lengths = [16, 17, 18, 19]
    indices = {}
    indices_found = False
    length_tuple_dict = collections.OrderedDict()
    for length in possible_lengths:
        length_tuple_dict[length] = []

    with open(data_file, 'r') as f:
        while not indices_found:
            line = f.readline().rstrip()
            if axes[0] and axes[1] in line:
                split_line = line.split()
                x_index = split_line.index(axes[0]) + 1
                y_index = split_line.index(axes[1]) + 1
                indices[axes[0]] = x_index
                indices[axes[1]] = y_index
                indices_found = True
        for line in f:
            split_line = line.rstrip().split()
            try:
                x = float(split_line[indices[axes[0]]])
                y = float(split_line[indices[axes[1]]])
                seq = split_line[0]
                seq_len = len(seq)
            except IndexError:
                continue
            if not math.isnan(x) and not math.isnan(y) and seq_len in possible_lengths:
                if seq_count_dict and seq in seq_count_dict.keys():
                    count = seq_count_dict[seq]
                    length_tuple_dict[seq_len].append((x, y, count))
                else:
                    length_tuple_dict[seq_len].append((x, y))
    return length_tuple_dict


def dataframe_from_file(data_file):
    columns = ['repressed', 'induced', 'free', 'fold_induction', 'fold_repression']
    with open(data_file, 'r') as f:
        dataframe = pd.read_table(f, sep='\s+', header=None, skiprows=2, index_col=0)
        dataframe.columns = columns
        dataframe.index.name = 'sequence'
    return dataframe


def add_column_for_coloring(dataframe, color_scheme):
    supported_schemes = ['length']
    if color_scheme == 'length':
        seq_list = list(dataframe.index)
        seq_length_list = []
        for x in seq_list:
            if type(x) is str:
                seq_length_list.append(len(x))
            else:
                seq_length_list.append(0)
        dataframe['length'] = pd.Series(seq_length_list, index=dataframe.index)
    else:
        raise Exception("color_scheme is not among one of the supported options. "
                        "the supported options are {0}".format(' '.join(supported_schemes)))
    return dataframe


def panda_plotly_plot(dataframe, axes, color_scheme, name, color_dict):
    repressed_foldinduction = [
            {
                'x': dataframe[dataframe['length'] == length].repressed,
                'y': dataframe[dataframe['length'] == length].fold_induction,
                'text': dataframe.index,
                'name': length,
                'mode': 'markers',
            } for length in [16, 17, 18, 19]
            ]
    repressed_foldinduction = go.Scatter(
        x=
    )
    free_foldinduction = [
        {
            'x': dataframe[dataframe['length'] == length].free,
            'y': dataframe[dataframe['length'] == length].fold_induction,
            'text': dataframe.index,
            'name': length,
            'mode': 'markers',
        } for length in [16, 17, 18, 19]
    ]

    update_menus = list([
        dict(active=-1,
             buttons=list([
                 dict(label='repression fold_induction',
                      method='update',
                      args=[{'visible': [True, False]},
                            {'title': 'CmeR'}]),

                 dict(label='free fold_induction',
                      method='update',
                      args=[{'visible': [False, True]},
                            {'title': 'CmeR'}])
             ]),
             )
    ])

    layout = dict(title='Test', showlegend=False, updatemenus=update_menus)
    print layout
    fig = {
        'data': [repressed_foldinduction, free_foldinduction],
        'layout': {
            'xaxis': {'title': axes[0], 'type': 'log'},
            'yaxis': {'title': axes[1]}
        }
    }
    url = plotly.offline.plot(fig, auto_open=False, show_link=False)


def plotly_data_files(data_file, x_axis, y_axis, color_scheme):
    plotly.plotly.sign_in('nwhoppe', 'ZKjPOgxRc3XFmeVrjlTJ')
    axes = [x_axis, y_axis]
    dataframe = dataframe_from_file(data_file)
    dataframe = add_column_for_coloring(dataframe, color_scheme)
    allowable_lengths = [16, 17, 18, 19]
    color_list = palette[len(allowable_lengths)]
    color_dict = dict(zip(allowable_lengths, color_list))
    dataframe = dataframe[dataframe[color_scheme].isin(allowable_lengths)]

    name = '{0}_{1}_{2}'.format(data_file.split('.')[0], x_axis, y_axis)
    # generate_plots(dataframe, axes, color_scheme, name, color_dict)
    panda_plotly_plot(dataframe, axes, color_scheme, name, color_dict)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of induction data from "
                                                 "genotype_fluorescence.py")
    parser.add_argument("-x", "--x_axis", default='repressed',
                        help="criterion to be plotted on x-axis (default: repressed)")
    parser.add_argument("-y", "--y_axis", default='fold_induction',
                        help="criterion to be plotted on y-axis (default: fold_induction)")
    parser.add_argument("-c", "--color_scheme", default='length',
                        help='specify how points on scatter plot should be colored. currently the supported options '
                             'are: length, ')

    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")
    args = parser.parse_args()
    plotly_data_files(args.data_file, args.x_axis, args.y_axis, args.color_scheme)
