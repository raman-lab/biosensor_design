#!/usr/bin/env python
import argparse
import pandas as pd
import plotly
import math
import numpy as np

from dna_shape_cluster import dataframe_from_file
from dnacurve import CurvedDNA
import plotly.graph_objs as go
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# result = CurvedDNA("ATCAGTGAACGCGCGTTCACTGAT", "TRINUCLEOTIDE", name="Example")
# res = result.curvature
#
# print np.size(res)
# print res
# print np.sum(res[0, :])


def scatter_plot_with_percentile(xy_dataframe, percentile, name):
    almost_gray = '#808080'
    almost_black = '#262626'
    # xy_dataframe = xy_dataframe[xy_dataframe['pssm score'] >= 0]
    x = np.asarray(xy_dataframe['curvature'])
    y = np.asarray(xy_dataframe['fold induction'])
    # x_bins = np.linspace(min(x), max(x), 16)
    # x_bins = range(int(min(x)) / 5, int(max(x)) + 5, 5)
    # bin_indices = np.digitize(x, x_bins)
    # y_percentile_medians = []
    #
    # for i in range(1, len(x_bins)):
    #     if i in bin_indices:
    #         y_bin = y[np.where(bin_indices == i)]
    #         q = np.percentile(y_bin, percentile)
    #         # y_above_percentile = [k for k in y_bin if k >= q]
    #         # y_above_percentile_median = np.median(y_above_percentile)
    #         # y_percentile_medians.append(y_above_percentile_median)
    #         y_percentile_medians.append(q)
    #     else:
    #         x_bins = np.delete(x_bins, i - 1)

    # x_avg_bins = []
    # for i in range(0, len(x_bins) - 1):
    #     avg = (x_bins[i] + x_bins[i + 1]) / 2.0
    #     x_avg_bins.append(avg)

    # x_smooth = np.linspace(min(x_avg_bins), max(x_avg_bins), 500)
    # y_smooth = scipy.interpolate.spline(x_avg_bins, y_percentile_medians, x_smooth)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15)
    # ax1.plot(x_smooth, y_smooth, color=almost_black, alpha=0.75)

    ax1.set_ybound(lower=0.0)
    # ax1.set_xbound(lower=min(x) - 1, upper=max(x) + 1)
    plt.xlabel('avg curvature angle', fontsize=20)
    plt.ylabel('fold induction', fontsize=20)

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

    plt.savefig('{0}.pdf'.format(name))
    plt.close()


def violin_plot(xy_dataframe, basename):
    pd.options.mode.chained_assignment = None
    # xy_dataframe = xy_dataframe[xy_dataframe['pssm score'] >= 0]
    # min_x = int(min(xy_dataframe['curve_angle'])) / 5
    # max_x = int(max(xy_dataframe['curve_angle'])) + 2
    min_x = 0.5 * (math.floor(min(xy_dataframe['curvature']) / 0.5))
    max_x = 0.5 * (math.ceil(max(xy_dataframe['curvature']) / 0.5))
    bins = np.arange(min_x, max_x + 0.5, 0.5)
    labels = []
    for i in range(0, len(bins) - 1):
        label = '{0}-{1}'.format(bins[i], bins[i + 1])
        labels.append(label)
    # print pd.cut(xy_dataframe.loc[:, 'pssm score'], bins=bins)
    xy_dataframe['bins'] = pd.cut(xy_dataframe['curvature'], bins=bins, labels=labels)
    sns.set_style("whitegrid")
    ax = sns.violinplot(x="bins", y="fold induction", data=xy_dataframe[['bins', 'fold induction']],
                        cut=0, scale='width', palette="muted")
    fig = ax.get_figure()
    fig.savefig('{0}_violin.pdf'.format(basename))


def dna_curvature_plotter(data_file, pssm_scores_json, allowable_lengths, column):
    pd.options.mode.chained_assignment = None
    dataframe = dataframe_from_file(data_file)
    # seq_list = list(dataframe.index)
    # seq_length_list = []
    # for x in seq_list:
    #     if type(x) is str:
    #         seq_length_list.append(len(x))
    #     else:
    #         seq_length_list.append(0)
    # dataframe['length'] = pd.Series(seq_length_list, index=dataframe.index)
    # dataframe = dataframe[dataframe['length'].isin(allowable_lengths)]
    # dataframe.dropna(subset=[column], how='any', inplace=True)
    # dataframe = dataframe[(dataframe.repressed < dataframe.repressed.quantile(0.8)) & (dataframe.fold_induction > 2)]

    # pssm_df = pd.read_json(pssm_scores_json)
    # pssm_df = pssm_df[pssm_df['pssm score'] >= 0]

    plotting_dict_3d = {
        16: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        },
        17: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        },
        18: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        },
        19: {
            'x': [],
            'y': [],
            'z': [],
            'text': []
        }
    }

    plotting_dict_2d = {
        16: {
            'x': [],
            'y': [],
            'text': []
        },
        17: {
            'x': [],
            'y': [],
            'text': []
        },
        18: {
            'x': [],
            'y': [],
            'text': []
        },
        19: {
            'x': [],
            'y': [],
            'text': []
        }
    }

    plotting_dict = {
        'fold induction': [],
        'curvature': [],
        'pssm score': [],
        'text': []
    }
    curvature_list = []
    length_data_dict = {16: [], 17: [], 18: [], 19: []}
    # for sequence in pssm_df.index:
    for sequence in dataframe.index:
        if not isinstance(sequence, float) and len(sequence) in length_data_dict.keys():
            seq_len = len(sequence)
            seq_curve_obj = CurvedDNA("TTGACA{0}TATAAT".format(sequence), "NUCLEOSOME", name="Example")
            curvature = seq_curve_obj.curvature[0, :]
            curvature = curvature[curvature.nonzero()]
            # curvature_sum = curvature.sum()
            # curvature_list.append(curvature_sum)
            # fold_ind = pssm_df.loc[sequence]['fold induction']
            # pssm_val = pssm_df.loc[sequence]['pssm score']
            # plotting_dict['curvature'].append(curvature_sum)
            # plotting_dict['pssm score'].append(pssm_val)
            # plotting_dict['fold induction'].append(fold_ind)
            # plotting_dict['text'].append(sequence)
            length_data_dict[seq_len].append(curvature)
        # plotting_dict_3d[seq_len]['x'].append(curvature)
        # plotting_dict_3d[seq_len]['y'].append(pssm_val)
        # plotting_dict_3d[seq_len]['z'].append(fold_ind)
        # plotting_dict_3d[seq_len]['text'].append('{0}'.format(sequence))
        #
        # plotting_dict_2d[seq_len]['x'].append(curvature)
        # plotting_dict_2d[seq_len]['y'].append(fold_ind)
        # plotting_dict_2d[seq_len]['text'].append('{0}'.format(sequence))

    # pssm_df['curvature'] = pd.Series(curvature_list, index=pssm_df.index)
    # new_df = pssm_df[['curvature', 'fold induction']].copy()
    # base_name = data_file.split('/')[-1].split('.')[0]
    # scatter_plot_with_percentile(new_df, 90, base_name)
    # violin_plot(new_df, base_name)

    base_name = '{0}'.format(data_file.split('/')[-1].split('.')[0])

    for length in length_data_dict:
        data = [
            go.Heatmap(
                # x=length_data_dict[length]['position'],
                # y=length_data_dict[length]['fold_induction'],
                z=length_data_dict[length],
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

    # data = []
    # trace = go.Scatter3d(
    #     x=plotting_dict['curvature'],
    #     y=plotting_dict['pssm score'],
    #     z=plotting_dict['fold induction'],
    #     text=plotting_dict['text'],
    #     mode='markers',
    #     opacity=0.8
    # )
    # data.append(trace)
    #
    # layout = go.Layout(
    #     title=base_name,
    #     scene=dict(
    #         xaxis=dict(
    #             title='curvature'
    #         ),
    #         yaxis=dict(
    #             title='pssm score'
    #         ),
    #         zaxis=dict(
    #             title='fold induction'
    #         ),
    #     ),
    # )
    #
    # fig = go.Figure(data=data, layout=layout)
    # plotly.offline.plot(fig, filename='{0}_curve_fi_pssm.html'.format(base_name), auto_open=False)
    #
    # fig = {
    #     'data': [
    #         {
    #             'x': plotting_dict['curvature'],
    #             'y': plotting_dict['fold induction'],
    #             'text': plotting_dict['text'],
    #             'mode': 'markers',
    #         }
    #         ],
    #     'layout': {
    #         'xaxis': {'title': 'curvature'},
    #         'yaxis': {'title': 'fold induction'},
    #         'title': base_name}
    # }
    #
    # plotly.offline.plot(fig, filename='{0}_curve_fi.html'.format(base_name), auto_open=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-l', '--lengths', nargs='*', default=[16, 17, 18, 19],
                        help='one heatmap will be made per length given')
    parser.add_argument('-col', '--column', default='fold_induction', help="column to use in data file")
    parser.add_argument('-p', '--pssm_scores_json',
                          help='json file with scores of sequences in data file to pssm profile of csi sequences')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-d', '--data_file', required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same as the x and y axes")

    args = parser.parse_args()
    dna_curvature_plotter(args.data_file, args.pssm_scores_json, args.lengths, args.column)