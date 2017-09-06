#!/usr/bin/env python

import argparse
import collections
import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from color_palettes import palette


def parse_plate_layout(plate_layout):
    """From plate layout file, create nested dictionaries - data - column header - each row - {OD,fluor,ratio}"""
    headers = collections.OrderedDict()
    with open(plate_layout, 'r') as f:
        lines = f.readline().split('\r')
        header_list = lines[0].split('\t')
    for header in header_list:
        headers[header] = collections.OrderedDict()
    for line in lines[1:]:
        split_line = line.split('\t')
        for c, column in enumerate(split_line):
            headers[headers.keys()[c]][column] = collections.OrderedDict([('od', []), ('fluor', [])])
    return headers


def is_time(time_stamp):
    """check formatting of time in first column"""
    try:
        time.strptime(time_stamp, '%H:%M:%S')
    except ValueError:
        return False
    else:
        return True


def convert_time_stamp(time_stamp):
    """round time to nearest half hour"""
    [hour, minute] = time_stamp.split(':')[:-1]
    nearest_half_hour = 0.5 * round(float(minute) / 30)
    return int(hour) + nearest_half_hour


def parse_plate_reader_data(data_file, data_dicts):
    """From data file, pull OD and fluor data to populate nested dictionaries. Assumes OD first then Fluor data"""
    times = collections.OrderedDict([('od', []), ('fluor', [])])
    with open(data_file, 'r') as f:

        lines = f.readline().split('\r')
    seen_first_time_data = False
    past_first_time_data = False
    for line in lines:
        split_line = line.split()
        try:
            time_stamp = split_line[0]
            row_data = split_line[2:]
        except IndexError:
            continue
        time_bool = is_time(time_stamp)

        if 'Time' in time_stamp and 'T' in split_line[1] and not seen_first_time_data:
            seen_first_time_data = True

        elif time_bool and not past_first_time_data and seen_first_time_data:
            times['od'].append(convert_time_stamp(time_stamp))

            jump_int = len(row_data) / len(data_dicts.keys())
            for h, header in enumerate(data_dicts.keys()):
                indices_to_take = range(h, len(row_data), jump_int)
                for c, column in enumerate(data_dicts[header].keys()):
                    data_to_add = row_data[indices_to_take[c]]
                    data_dicts[header][column]['od'].append(data_to_add)

        elif 'Time' in time_stamp and 'T' in split_line[1] and seen_first_time_data:
            past_first_time_data = True

        elif time_bool and past_first_time_data and len(times['fluor']) < len(times['od']):
            times['fluor'].append(convert_time_stamp(time_stamp))

            jump_int = len(row_data) / len(data_dicts.keys())
            for h, header in enumerate(data_dicts.keys()):
                indices_to_take = range(h, len(row_data), jump_int)
                for c, column in enumerate(data_dicts[header].keys()):
                    data_to_add = row_data[indices_to_take[c]]
                    data_dicts[header][column]['fluor'].append(data_to_add)

        else:
            continue
    return times, data_dicts


def plot_plate_reader_data(time_data, data_dicts):
    """makes png scatter plot of data: od, fluor, and ratio of fluor:od"""
    almost_gray = '#808080'
    almost_black = '#262626'
    yaxis_dict = {'od': 'OD600', 'fluor': 'Fluorescence (RFU)', 'ratio': 'Fluorescence / OD600'}
    color_set = [almost_gray]
    color_set.extend(palette[len(data_dicts[data_dicts.keys()[0]].keys()) - 1])
    time_data = time_data['fluor']
    for d_type in ['od', 'fluor', 'ratio']:
        for header in data_dicts.keys():
            fig = plt.figure()

            for c, column in enumerate(data_dicts[header].keys()):
                if d_type is 'ratio':
                    od_data = np.array(data_dicts[header][column]['od'], dtype=float)
                    fluor_data = np.array(data_dicts[header][column]['fluor'], dtype=float)
                    dtype_data = np.divide(fluor_data, od_data)
                else:
                    dtype_data = np.array(data_dicts[header][column][d_type], dtype=float)

                ax1 = fig.add_subplot(111)
                ax1.scatter(time_data, dtype_data, c=color_set[c], marker='o', alpha=0.9, edgecolor=almost_black,
                            linewidth=0.15, label='{0}'.format(column))

            legend = ax1.legend(loc=2, scatterpoints=1, framealpha=1)
            rect = legend.get_frame()
            rect.set_linewidth(0.25)
            texts = legend.texts
            for t in texts:
                t.set_color(almost_black)
            plt.xlim(xmin=0, xmax=max(time_data) + 1)
            plt.xlabel('Time (hours)', fontsize=20)
            plt.ylabel(yaxis_dict[d_type], fontsize=20)

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
            plt.savefig('{0}_{1}.png'.format(header, d_type), format='png', dpi=1000)
            plt.close()


def plate_reader_data(data, plate_layout):
    plate_data_dicts = parse_plate_layout(plate_layout)
    time_data, data_dicts = parse_plate_reader_data(data, plate_data_dicts)
    plot_plate_reader_data(time_data, data_dicts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="script to make scatter plots of OD, fluorescence, and fluorescence/OD"
                                                 " data from plate reader. script is specific to output formatting of "
                                                 "Raman Lab's plate reader. script requires plain text "
                                                 "plate reader data as well as a plain text plate "
                                                 "layout file as input (export option in excel).")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-d", "--data", required=True,
                           help="one plain text data file from plate reader")
    requiredO.add_argument("-p", "--plate_layout", required=True,
                           help='one plain text file. columns correspond to different plots. first row is '
                                'header. subsequent rows used to label data on plots')

    args = parser.parse_args()
    plate_reader_data(args.data, args.plate_layout)
