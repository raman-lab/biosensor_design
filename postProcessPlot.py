#!/usr/bin/env python

# made for comparing unfiltered and filtered scorefiles for Rosetta enzdes post analysis

import argparse
import collections
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages       


def data_from_sc_file(axes, f, uf, true_max):
    """initializes two dictionaries and poplulates them based on -f and -u options"""
    f_combo_dict = collections.defaultdict(list)
    uf_combo_dict = collections.defaultdict(list)
    max_x = -10000
    max_y = -10000
    min_x = 10000
    min_y = 10000

    for fileType in [uf, f]:
        for i, item in enumerate(fileType):
            with open(item) as f:
                header = f.readline().split()
                indices = [header.index(a) for a in axes]
                
                for line in f:
                    line_list = line.split()
                    if (not line_list) or (line_list[0].startswith("#")) or (line_list[0][0].isalpha()):
                        continue
                    try:
                        desc_str = line_list[indices[-1]]
                        found_desc = re.search('A([0-9]+)_P([0-9]+)', desc_str).group()
                    except AttributeError:
                        continue
                    
                    point_list = [line_list[i] for i in indices[:-1]]
                    point_tuple = tuple(map(float, point_list))
                    if point_tuple[0] > max_x:
                        max_x = point_tuple[0]
                    if point_tuple[0] < min_x:
                        min_x = point_tuple[0]
                    if point_tuple[1] > max_y:
                        max_y = point_tuple[1]
                    if point_tuple[1] < min_y:
                        min_y = point_tuple[1]

                    if not true_max:
                        if max_x > 0:
                            max_x = 0
                        if max_y > 0:
                            max_y = 0
                    
                    if fileType == uf:
                        uf_combo_dict[found_desc].append(point_tuple)
                    else:
                        f_combo_dict[found_desc].append(point_tuple)
    return uf_combo_dict, f_combo_dict, min_x, max_x, min_y, max_y


def gen_plots(uf_dict, f_dict, min_x, max_x, min_y, max_y, axes, name, histogram, total):
    """makes pdf of plots - one plot for each A[0-9]_P[0-9]"""
    with PdfPages(name) as pdf:
        total_xuf = []
        total_yuf = []
        total_xf = []
        total_yf = []
        for entry in uf_dict:
            print 'Making plot for ' + entry
            xuf, yuf = zip(*uf_dict[entry])
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.scatter(xuf, yuf, c='#ad4851', marker='o', label='initial structures')
            try:
                xf, yf = zip(*f_dict[entry])
                ax1.scatter(xf, yf, c='orange', marker='x', label='selected structures')
            except ValueError:
                xf = []
                yf = []
            plt.legend(loc='upper right')
            plt.title(entry, fontsize=30)
            plt.xlim(min_x, max_x)
            plt.ylim(min_y, max_y)
            plt.xlabel(axes[0], fontsize=20)
            plt.ylabel(axes[1], fontsize=20)
            pdf.savefig(fig)
            plt.close()

            if total:
                total_xuf.extend(xuf)
                total_yuf.extend(yuf)
                total_xf.extend(xf)
                total_yf.extend(yf)

            if histogram:
                bins = np.linspace(min_y, max_y, num=10)
                plt.hist(yuf, bins, alpha=0.5, color='b', label='initial structures')
                try:
                    plt.hist(yf, bins, alpha=0.5, color='orange', label='selected structures')
                except ValueError:
                    pass
                plt.legend(loc='upper right')
                plt.title(entry, fontsize=30)
                plt.xlabel(axes[1], fontsize=20)
                plt.ylabel('Frequency', fontsize=20)
                pdf.savefig()
                plt.close()

        if total:
            print 'Making composite plot'
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.scatter(total_xuf, total_yuf, c='#ad4851', marker='o', label='initial structures')
            ax1.scatter(total_xf, total_yf, c='orange', marker='x', label='selected structures')
            plt.legend(loc='upper right')
            plt.title('Composite Plot', fontsize=30)
            plt.xlim(min_x, max_x)
            plt.ylim(min_y, max_y)
            plt.xlabel(axes[0], fontsize=20)
            plt.ylabel(axes[1], fontsize=20)
            pdf.savefig(fig)
            plt.close()


def main(x_axis, y_axis, filtered, unfiltered, name, histogram, total, true_max):
    """create axes variable and calls previous functions"""
    axes = [x_axis, y_axis, 'description']
    uf_dict, f_dict, min_x, max_x, min_y, max_y = data_from_sc_file(axes, filtered, unfiltered, true_max)
    gen_plots(uf_dict, f_dict, min_x, max_x, min_y, max_y, axes, name, histogram, total)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generates scatter plot of data from rosetta score files")
    parser.add_argument("-x", "--xaxis",
                        help="criterion to be plotted on x-axis (default: total_score)",
                        default='total_score')
    parser.add_argument("-y", "--yaxis",
                        help="criterion to be plotted on y-axis (default: SR_1_total_score)",
                        default='SR_1_total_score')
    parser.add_argument("-n", "--name", default='postProcessPlot.pdf',
                        help='name of output pdf (default: postProcessPlot.pdf')
    parser.add_argument("-b", "--histogram", action="store_true",
                        help="turn on histogram for y-axis parameter")
    parser.add_argument("-c", "--composite", action="store_true",
                        help='make a composite plot that combines all subplots')
    parser.add_argument("-t", "--true_max", action="store_true",
                        help='make plots with true maximum - will not cap max at 0')
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-s", "--selected", nargs='*', required=True,
                           help="one or more filtered score files from which data is pulled")
    requiredO.add_argument("-i", "--initial", nargs='*', required=True,
                           help="one or more unfiltered score files from which data is pulled")
    args = parser.parse_args()

    main(args.xaxis, args.yaxis, args.selected, args.initial, args.name, args.histogram, args.composite, args.true_max)

    


