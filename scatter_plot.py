#!/usr/bin/env python
import argparse
import collections
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
from color_palettes import palette


def get_scores(axes, score_file):
    """get data specified by axes options from the score file and store as tuples in a list"""
    data_list = []
    with open(score_file, 'r') as f:
        header = f.readline().split()
        indices = [header.index(a) for a in axes]

        for line in f:
            line_list = line.split()
            if (not line_list) or (line_list[0].startswith("#")) or (line_list[0][0].isalpha()):
                continue

            point_list = [line_list[i] for i in indices[:-1]]
            point_tuple = tuple(map(float, point_list))
            data_list.append(point_tuple)
    return data_list


def generate_plots(data_tuple_dict, axes, name, legend_list):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[len(data_tuple_dict) - 1]
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

    x_min = float(sys.maxint)
    y_min = float(sys.maxint)
    with PdfPages(name) as pdf:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        for i, score_file in enumerate(data_tuple_dict.keys()):

            x, y = zip(*data_tuple_dict[score_file])
            if i is 0:
                ax1.scatter(x, y, c=almost_gray, marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                            label='{0}'.format(legend_list[i]))
            else:
                ax1.scatter(x, y, c=color_set[i-1], marker='o', alpha=0.5, edgecolor=almost_black, linewidth=0.15,
                            label='{0}'.format(legend_list[i]))

            if min(x) < x_min:
                x_min = min(x)
            if min(y) < y_min:
                y_min = min(y)

        legend = ax1.legend(loc='best', scatterpoints=1, framealpha=0.5)
        rect = legend.get_frame()
        rect.set_linewidth(0.0)
        texts = legend.texts
        for t in texts:
            t.set_color(almost_black)
        plt.xlim(x_min, xmax=0)
        plt.ylim(y_min, ymax=0)
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


def plot_score_files(x_axis, y_axis, score_files, name, description_list):
    """create axes variable and calls previous functions"""
    axes = [x_axis, y_axis, 'description']
    data_tuple_lists = collections.OrderedDict()
    for score_file in score_files:
        data_tuple_lists[score_file] = get_scores(axes, score_file)
    generate_plots(data_tuple_lists, axes, name, description_list)


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

    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-s", "--score_files", nargs='*', required=True,
                           help="one or more filtered score files from which data is pulled")
    requiredO.add_argument("-d", "--description", nargs='*', required=True,
                           help='provide short description of each score file in the same order that the score '
                                'files are given. the description will be used in the legend of the plot.')

    args = parser.parse_args()

    plot_score_files(
        args.x_axis, args.y_axis, args.score_files, args.name, args.description
    )
