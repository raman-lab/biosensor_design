#!/usr/bin/env python
import argparse
import collections
import itertools
from color_palettes import palette
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def generate_line_plots(file_len_data_dict, name):
    almost_gray = '#808080'
    almost_black = '#262626'
    with PdfPages('{0}.pdf'.format(name)) as pdf:
        for dna_shape_file, tag_data_dict in file_len_data_dict.iteritems():
            pro = dna_shape_file.split('.')[0]
            metric = dna_shape_file.split('.')[-1]
            color_set = palette[len(tag_data_dict)]

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            for i, tag in enumerate(tag_data_dict.keys()):
                y_list = tag_data_dict[tag]
                ax1.plot(range(1, len(y_list) + 1), y_list, '-o', c=color_set[i],
                         alpha=0.5, linewidth=0.5, label='{0}'.format(tag.split('>')[-1].split()[-1]))

            legend = ax1.legend(bbox_to_anchor=(1.04, 0.5), loc='center left', scatterpoints=1, framealpha=1)
            rect = legend.get_frame()
            rect.set_linewidth(0.25)
            texts = legend.texts
            for t in texts:
                t.set_color(almost_black)
            plt.ylabel('{0}'.format(metric), fontsize=20)
            plt.xlabel('position', fontsize=20)
            plt.title('{0}_{1}'.format(pro, metric))
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
            pdf.savefig(fig, bbox_extra_artists=(legend,), bbox_inches='tight')
            plt.close()


def dna_shape_line_plotter(dna_shape_input_files):
    # tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta_file)
    file_tag_data_dict = collections.defaultdict(dict)
    for dna_shape_file in dna_shape_input_files:
        name = dna_shape_file.split('.')[0]
        file_tag_data_dict[dna_shape_file] = {}
        with open(dna_shape_file, 'r') as f:
            for tag, data_string in itertools.izip_longest(f, f, fillvalue=None):
                tag = tag.rstrip().split('>', 1)[-1]
                string_list = [x for x in data_string.rstrip().split(',') if x != 'NA']
                data_list = map(float, string_list)
                # seq = tag_seq_dict[tag]
                file_tag_data_dict[dna_shape_file][tag] = data_list
    generate_line_plots(file_tag_data_dict, name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of dna shape parameters")
    required = parser.add_argument_group('required arguments')
    required.add_argument('-s', '--dna_shape_files', nargs='*', required=True,
                          help='files from DNAshape. at least one file must be provided')
    parser.add_argument('-n', '--name')
    args = parser.parse_args()
    dna_shape_line_plotter(args.dna_shape_files)
