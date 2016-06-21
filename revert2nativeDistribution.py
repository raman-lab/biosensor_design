#!/usr/bin/env python


import argparse
import numpy as np
import matplotlib.pyplot as plt


def main(input_files, histogram, plot_name):
    res_before_lst = []
    res_after_lst = []
    for in_file in input_files:
        with open(in_file) as f:
            for line in f:
                if line.startswith('Sequences'):
                    line_list = line.split()
                    res_before = float(line_list[4])
                    res_after = float(line_list[9])
                    res_before_lst.append(res_before)
                    res_after_lst.append(res_after)
    avg_before = np.around(np.average(res_before_lst), decimals=2)
    var_before = np.around(np.std(res_before_lst), decimals=2)
    avg_after = np.around(np.average(res_after_lst), decimals=2)
    var_after = np.around(np.std(res_after_lst), decimals=2)
    print "The average number of designed amino acids per structure is {0} with stdev {1} before reversion".format(
        avg_before, var_before)
    print "The average number of designed amino acids per structure is {0} with stdev {1} after reversion".format(
        avg_after, var_after)
    if histogram:
        bins = np.linspace(0, 15, 16)
        plt.hist(res_before_lst, bins, alpha=0.5, label='{0} +/- {1} mutations before'.format(
            avg_before, var_before))
        plt.hist(res_after_lst, bins, alpha=0.5, label='{0} +/- {1} mutations after'.format(
            avg_after, var_after))
        plt.legend(loc='upper right')
        if not plot_name:
            plot_name = 'revert_to_native_distribution.png'
        plt.title(plot_name)
        plt.savefig(plot_name, format='png', dpi=1000)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Script to calculate distribution of reverted amino acids and display histogram, if specified,
        given log files from the revert_to_native application in Rosetta""")
    parser.add_argument("-p", "--plot_histogram", action="store_true",
                        help="turn on histogram")
    parser.add_argument("-n", "--name",
                        help="name of output file with .png extension (default: revert_to_native_distribution.png")
    parser.add_argument("log_files", nargs='*', help='one or more log files from revert_to_native run')
    args = parser.parse_args()
    main(args.log_files, args.plot_histogram, args.name)
