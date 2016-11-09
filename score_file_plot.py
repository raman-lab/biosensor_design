#!/usr/bin/env python
import argparse
import collections
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from color_palettes import palette


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
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(uf_dict)]
    common_score_terms = {
        'total_score': 'Total Energy (REU)',
        'fa_rep': 'Total Repulsive Energy (REU)',
        'hbond_sc': 'Total Hydrogen Bond Side Chain Energy (REU)',
        'all_cst': 'Total Constraint Energy (REU)',
        'tot_pstat_pm': 'Pack Statistics with Ligand',
        'tot_nlpstat_pm': 'Pack Statistics without Ligand',
        'tot_burunsat_pm': 'Total Number of Buried Unsatisfied Polar Residues',
        'tot_hbond_pm': 'Total Number of Hydrogen Bonds',
        'tot_NLconts_pm': 'Total Number of Non-local Contacts',
        'SR_1_total_score': 'Total Ligand Energy (REU)',
        'SR_1_fa_rep': 'Ligand Repulsive Energy (REU)',
        'SR_1_hbond_sc': 'Ligand Hydrogen Bond Energy (REU)',
        'SR_1_all_cst': 'Ligand Constraint Energy (REU)',
        'SR_1_hbond_pm': 'Number of Ligand Hydrogen Bonds (REU)',
        'SR_1_burunsat_pm': 'Number of Buried Unsatisfied Polar Residues in Binding Pocket'
    }

    with PdfPages(name) as pdf:
        total_xuf = []
        total_yuf = []
        total_xf = []
        total_yf = []
        for index, entry in enumerate(uf_dict):
            print 'Making plot for ' + entry
            xuf, yuf = zip(*uf_dict[entry])
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.scatter(xuf, yuf, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                        label='before reversion')
            try:
                xf, yf = zip(*f_dict[entry])
                ax1.scatter(xf, yf, c=color_set[index], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                            label='after reversion')
            except ValueError:
                xf = []
                yf = []

            legend = ax1.legend(loc='best', scatterpoints=1, framealpha=0.5)
            rect = legend.get_frame()
            rect.set_linewidth(0.0)
            texts = legend.texts
            for t in texts:
                t.set_color(almost_black)
            plt.title('Distribution of Reverted {0} Energies'.format(entry), fontsize=25, y=1.02)
            plt.xlim(min_x, max_x)
            plt.ylim(min_y, max_y)
            plt.xlabel(common_score_terms[axes[0]], fontsize=20)
            plt.ylabel(common_score_terms[axes[1]], fontsize=20)

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

            if total:
                total_xuf.extend(xuf)
                total_yuf.extend(yuf)
                total_xf.extend(xf)
                total_yf.extend(yf)

            if histogram:
                bins = np.linspace(min_y, max_y, num=10)
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.hist(yuf, bins, alpha=0.5, color=almost_gray, label='before reversion')
                try:
                    plt.hist(yf, bins, alpha=0.5, color=color_set[index], label='after reversion')
                except ValueError:
                    pass

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

                plt.title('Distribution of Reverted {0} Energies'.format(entry), fontsize=25, y=1.02)
                plt.xlabel(common_score_terms[axes[1]], fontsize=20)
                plt.ylabel('Frequency', fontsize=20)
                pdf.savefig()
                plt.close()

        if total:
            print 'Making composite plot'
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.scatter(total_xuf, total_yuf, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black,
                        linewidth=0.15, label='all structures before filter')
            for index, entry in enumerate(uf_dict):
                try:
                    xf, yf = zip(*f_dict[entry])
                    ax1.scatter(xf, yf, c=color_set[index], marker='o', alpha=0.5, edgecolor=almost_black,
                                linewidth=0.15, label='all structures after filter')
                except ValueError:
                    pass

            # legend = ax1.legend(loc='best', scatterpoints=1, framealpha=0.5)
            # rect = legend.get_frame()
            # rect.set_linewidth(0.0)
            # texts = legend.texts
            # for t in texts:
            #     t.set_color(almost_black)

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

            plt.title('Distribution of Reverted Energies', fontsize=25, y=1.02)
            plt.xlim(min_x, max_x)
            plt.ylim(min_y, max_y)
            plt.xlabel(common_score_terms[axes[0]], color=almost_black, fontsize=20)
            plt.ylabel(common_score_terms[axes[1]], color=almost_black, fontsize=20)
            pdf.savefig(fig)
            plt.close()


def plot_scores(x_axis, y_axis, selected, initial, name, histogram, composite, true_max):
    """create axes variable and calls previous functions"""
    axes = [x_axis, y_axis, 'description']
    uf_dict, f_dict, min_x, max_x, min_y, max_y = data_from_sc_file(axes, selected, initial, true_max)
    gen_plots(uf_dict, f_dict, min_x, max_x, min_y, max_y, axes, name, histogram, composite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of data from rosetta score files")
    parser.add_argument("-x", "--x_axis",
                        help="criterion to be plotted on x-axis (default: total_score)",
                        default='total_score')
    parser.add_argument("-y", "--y_axis",
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

    plot_scores(
        args.x_axis, args.y_axis, args.selected, args.initial,
        args.name, args.histogram, args.composite, args.true_max
    )
