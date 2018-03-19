#!/usr/bin/env python
import argparse
import collections
import itertools
import json
from color_palettes import palette
import math
import matplotlib
from induction_cluster_plotly import parse_fasta_to_tag_seq_dict
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import scipy.optimize


def fit_sin(tt, yy):
    """Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq",
    period" and 'fitfunc'"""
    tt = numpy.array(tt)
    yy = numpy.array(yy)
    ff = numpy.fft.rfftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(numpy.fft.rfft(yy))
    guess_freq = abs(ff[numpy.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = numpy.std(yy) * 2.**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.*numpy.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c): return A * numpy.sin(w*t + p) + c
    try:
        popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    except RuntimeError:
        return False
    A, w, p, c = popt
    f = w/(2.*numpy.pi)
    fitfunc = lambda t: A * numpy.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": numpy.max(pcov), "rawres": (guess,popt,pcov)}


def generate_plots(file_len_data_dict):
    """makes pdf scatter plot"""
    # set up colors and names
    almost_gray = '#808080'
    almost_black = '#262626'
    color_set = palette[4]
    log_axes = ['repressed', 'induced', 'free']
    normal_axes = ['fold_induction', 'fold_repression', 'csi_score']

    for dna_shape_file, length_data_dict in file_len_data_dict.iteritems():
        pro = dna_shape_file.split('.')[0]
        metric = dna_shape_file.split('.')[-1]
        name = '{0}_{1}.pdf'.format(pro, metric)
        with PdfPages(name) as pdf:
            for param in ['amp', 'omega', 'offset']:
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                cumulative_list = []
                for i, length in enumerate(length_data_dict.keys()):
                    sin_dict_list, y_list = zip(*length_data_dict[length])
                    param_list = []
                    for sin_dict in sin_dict_list:
                        param_list.append(sin_dict[param])
                    cumulative_list.extend(param_list)
                    ax1.scatter(param_list, y_list, marker='o', c=color_set[i-1],
                                alpha=0.5, edgecolor=almost_black, linewidth=0.15, label='{0}'.format(length))
                values = numpy.array(cumulative_list)
                med = numpy.median(values)
                mad = numpy.median(numpy.abs(values - med))
                x_lower = med - 10 * mad
                x_upper = med + 10 * mad

                if x_lower < min(cumulative_list):
                    x_lower = min(cumulative_list)
                if x_upper > max(cumulative_list):
                    x_upper = max(cumulative_list)

                ax1.set_xbound(lower=x_lower, upper=x_upper)
                legend = ax1.legend(loc='best', scatterpoints=1, framealpha=1)
                rect = legend.get_frame()
                rect.set_linewidth(0.25)
                texts = legend.texts
                for t in texts:
                    t.set_color(almost_black)
                plt.ylabel('fold induction', fontsize=20)
                plt.xlabel('{0}'.format(param), fontsize=20)
                plt.title('{0}_{1}_{2}'.format(pro, metric, param))
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


def get_data_by_column(data_file, axis):
    possible_lengths = [16, 17, 18, 19]
    indices = {}
    indices_found = False
    length_tuple_dict = collections.OrderedDict()
    for length in possible_lengths:
        length_tuple_dict[length] = {}
    with open(data_file, 'r') as f:
        while not indices_found:
            line = f.readline().rstrip()
            if axis in line:
                split_line = line.split()
                index = split_line.index(axis) + 1
                indices[axis] = index
                indices_found = True
        for line in f:
            split_line = line.rstrip().split()
            try:
                y = float(split_line[indices[axis]])
                seq = split_line[0]
                seq_len = len(seq)
            except IndexError:
                continue
            if not math.isnan(y) and seq_len in possible_lengths:
                length_tuple_dict[seq_len][seq] = y
    return length_tuple_dict


def dna_shape_plotter(data_file, y_axis, dna_shape_input_files, fasta_file):
    len_seq_data_dict = get_data_by_column(data_file, y_axis)
    tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta_file)
    file_len_data_dict = collections.defaultdict(dict)
    for dna_shape_file in dna_shape_input_files:
        file_len_data_dict[dna_shape_file] = {16: [], 17: [], 18: [], 19: []}
        with open(dna_shape_file, 'r') as f:
            for tag, data_string in itertools.izip_longest(f, f, fillvalue=None):
                tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
                # string_list = ['nan' if x == 'NA' else x for x in data_string.rstrip().split(',')]
                string_list = [x for x in data_string.rstrip().split(',') if x != 'NA']
                data_list = map(float, string_list)
                seq = tag_seq_dict[tag]
                length = len(seq)
                if length in file_len_data_dict[dna_shape_file].keys():
                    y_data = len_seq_data_dict[length][seq]
                    sin_dict = fit_sin(range(len(data_list)), data_list)
                    if sin_dict:
                        data_tup = (sin_dict, y_data)
                        file_len_data_dict[dna_shape_file][length].append(data_tup)
    generate_plots(file_len_data_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script generates scatter plot of dna shape parameters")
    required = parser.add_argument_group('required arguments')
    required.add_argument("-d", "--data_file", required=True,
                          help="columns of data separated by white space. columns used to generate plot must be "
                               "named the same color_axis")
    required.add_argument('-y', '--y_axis', default='fold_induction',
                          help='data for y axis. default is "fold_induction"')
    required.add_argument('-s', '--dna_shape_files', nargs='*', required=True,
                          help='files from DNAshape. at least one file must be provided')
    required.add_argument('-f', '--fasta', required=True,
                          help='fasta file from data table that was used as input to dna shape')
    args = parser.parse_args()
    dna_shape_plotter(args.data_file, args.y_axis, args.dna_shape_files, args.fasta)
