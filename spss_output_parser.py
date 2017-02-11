#!/usr/bin/env python

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt


def spss_parse_data(spss_file):
    descriptive_table = ''
    post_hoc_tests = ''

    with open(spss_file, 'rU') as f:
        data_chunks = f.read().split(os.linesep + os.linesep)

    for i, chunk in enumerate(data_chunks):
        if chunk.startswith('Descriptives'):
            descriptive_table = chunk
        elif 'Post Hoc Tests' in chunk:
            post_hoc_tests = chunk

    if (not descriptive_table) and (not post_hoc_tests):
        raise Exception("input file may not be formatted as expected. attempted to use 'Descriptives' and "
                        "'Post Hoc Tests' to grab their associated tables, but got an empty string")

    return descriptive_table, post_hoc_tests


def break_apart_table(spss_table):
    line_list = []
    for line in spss_table.splitlines():
        if line.startswith('|') and ('--' not in line):
            split_line = line.split('|')
            new_line = []
            for piece in split_line:
                new_line.append(piece.strip())
            line_list.append(new_line)
    return line_list


def spss_parse_descriptive(descriptive_table):
    labels, mean, stdev = [[], [], []]
    line_list = break_apart_table(descriptive_table)

    try:
        mean_index = line_list[0].index('Mean')
    except ValueError:
        raise ValueError("Could not find 'Mean' column in descriptives table")
    try:
        stdev_index = line_list[0].index('Std. Deviation')
    except ValueError:
        raise ValueError("Could not find 'Std. Deviation' in descriptives table")

    for members in line_list[2:-1]:
        if members[1] is not 'Total':
            labels.append(members[1])
            mean.append(float(members[mean_index]))
            stdev.append(float(members[stdev_index]))

    if not labels:
        raise Exception("descriptives table may not be formatted as expected. could not find labels in first column")
    if not mean:
        raise Exception("descriptives table may not be formatted as expected. "
                        "could not find means in column labeled 'Means'")
    if not stdev:
        raise Exception("descriptives table may not be formatted as expected. "
                        "could not find stdev in column labeled 'Std. Deviation'")
    if not len(labels) == len(mean) == len(stdev):
        raise Exception("descriptives table may not be formatted as expected. "
                        "number of labels, means, and stdevs do not match")
    return labels, mean, stdev


def spss_parse_post_hoc(post_hoc_table, first_descriptive_label):
    significant = []
    line_list = break_apart_table(post_hoc_table)
    try:
        mean_diff_index = line_list[0].index('Mean Difference (I-J)')
    except ValueError:
        raise ValueError("Could not find 'Mean Difference (I-J)' column in Post Hoc table")
    past_first_descriptive_label = False
    for members in line_list[2:]:
        if (members[1] != first_descriptive_label) and (members[1]):
            past_first_descriptive_label = True
        if (not past_first_descriptive_label) and ('*' in members[mean_diff_index]):
            significant.append(members[2])
    return significant


def spss_bar_chart(labels, mean, stdev, significant_entries, filename):
    almost_gray = '#808080'
    almost_black = '#262626'
    fig, ax1 = plt.subplots()

    index = np.arange(len(labels))
    bar_width = 0.5

    rects1 = plt.bar(index, mean, bar_width,
                     alpha=0.5,
                     color=almost_gray,
                     yerr=stdev,
                     ecolor=almost_black,
                     capsize=5
                     )

    asterisks = [None] * len(labels)
    for entry in significant_entries:
        entry_index = labels.index(entry)
        asterisks[entry_index] = '*'

    for i, asterisk in enumerate(asterisks):
        if asterisk is None:
            continue
        else:
            ax1.text(i + 0.1, mean[i] + .1, asterisk, color=almost_black, fontweight='bold')

    plt.ylabel('Mean')
    plt.xticks(index + bar_width / 2, labels, fontsize='8')
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

    plt.savefig("{0}.png".format(filename.split('.')[0]), format='png', dpi=1000)
    plt.close()


def spss_parse(spss_file_list):
    for spss_file in spss_file_list:
        descriptive_table, post_hoc_table = spss_parse_data(spss_file)
        labels, mean, stdev = spss_parse_descriptive(descriptive_table)
        significant_entries = spss_parse_post_hoc(post_hoc_table, labels[0])
        spss_bar_chart(labels, mean, stdev, significant_entries, spss_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to parse the output files from spss software and create bar charts from data.
        grabs data from descriptive table and significance from first sub-table in Post Hoc Tests table"""
    )
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-f", "--file", nargs='*', required=True,
                           help="one or more spss text files")

    args = parser.parse_args()
    spss_parse(args.file)
