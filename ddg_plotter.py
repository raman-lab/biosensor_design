#!/usr/bin/env python

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def ddg_histogram(ddg_scores, ddg_mean, ddg_stdev, plot_name):
    almost_gray = '#808080'
    almost_black = '#262626'
    fig, ax1 = plt.subplots()
    max_bin = max(ddg_scores)
    if max_bin > (ddg_mean + 3 * ddg_stdev):
        max_bin = ddg_mean + 3 * ddg_stdev
    bins = np.linspace(0, max_bin, num=20)
    ax1.hist(ddg_scores, bins, alpha=0.5, color=almost_gray, label='{0} +/- {1}'.format(ddg_mean, ddg_stdev))

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

    plt.savefig("{0}.png".format(plot_name), format='png', dpi=1000)
    plt.close()


def is_number(input_string):
        try:
            float(input_string)
        except ValueError:
            return False
        else:
            return True


def ddg_get_scores(out_file):
    scores = []
    with open(out_file, 'r') as f:
        for line in f:
            line_list = line.split()
            if len(line) < 2:
                continue
            elif is_number(line_list[2]):
                scores.append(float(line_list[2]))
    return scores


def ddg_get_stats(ddg_scores):
    return np.mean(ddg_scores), np.std(ddg_scores)


def ddg_plotter(out_files):
    for out_file in out_files:
        ddg_scores = ddg_get_scores(out_file)
        ddg_mean, ddg_stdev = ddg_get_stats(ddg_scores)
        plot_name = out_file.split('.')[0]
        ddg_histogram(ddg_scores, ddg_mean, ddg_stdev, plot_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""script to make histogram and calculate mean and standard deviation from Rosetta ddg
        calculation output files"""
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-o', '--out_files', nargs='*', required=True,
                          help='output file(s) from Rosetta ddg calculations.')
    args = parser.parse_args()
    ddg_plotter(args.out_files)
