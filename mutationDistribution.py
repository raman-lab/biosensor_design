#!/usr/bin/env python

# Calculates number of mutations per structure and makes heat map
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def compare_seq(wt_fasta, query_fasta_list):
    """populates dictionary of mutated positions
    key is a tuple (residue, position) and key is frequency count"""
    mut_pos_dict = {}
    num_mut_per_query = []
    with open(wt_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            else:
                wt_seq = line

    for query in query_fasta_list:
        with open(query, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                else:
                    query_seq = line
                mut_count = 0
                for r, residue in enumerate(query_seq):
                    if residue is wt_seq[r]:
                        continue
                    else:
                        mut_count += 1
                        try:
                            mut_pos_dict[(r+1, residue)] += 1
                        except KeyError:
                            mut_pos_dict[(r+1, residue)] = 1

                num_mut_per_query.append(mut_count)
    return mut_pos_dict, num_mut_per_query, wt_seq


def dict2array3d(d):
    """extract numpy arrays from dictionary"""
    sorted_list = sorted(d.items())
    key_tup, value_tup = zip(*sorted_list)
    key0_tup, key1_tup = zip(*key_tup)
    value_array = np.array(value_tup)
    key0_array = np.array(key0_tup)
    key1_array = np.array(key1_tup)
    return key0_array, key1_array, value_array


def mk_heat_map(mut_dict, wt_seq, query_fasta_list, pwm_bool, resfile, name):
    """create heat map of mutated residues and use subplot to identify wt residues"""
    pos_array, res_array, count_array = dict2array3d(mut_dict)
    almost_black = '#262626'

    unique_mut_pos_list = list(np.unique(pos_array))
    wt_dict = {}
    for query in query_fasta_list:
        with open(query, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    continue
                else:
                    query_seq = line
                    for p in unique_mut_pos_list:
                        if query_seq[p-1] is wt_seq[p-1]:
                            try:
                                wt_dict[(p, wt_seq[p-1])] += 1
                            except KeyError:
                                wt_dict[(p, wt_seq[p-1])] = 1
    mut_dict.update(wt_dict)
    tot_pos_array, tot_res_array, tot_count_array = dict2array3d(mut_dict)

    unique_pos_array = np.unique(tot_pos_array)
    unique_res_array = np.unique(tot_res_array)
    xi = np.arange(len(unique_pos_array))
    yi = np.arange(len(unique_res_array))
    Xi, Yi = np.meshgrid(xi, yi)

    Z = np.zeros(Xi.shape)
    xbox = []
    ybox = []
    if pwm_bool:
        cmap_color = 'coolwarm'
        for i in xi:
            indices = np.where(tot_pos_array == unique_pos_array[i])[0]
            N = sum(tot_count_array[indices[:]])
            for j in yi:
                try:
                    p = float(mut_dict[(unique_pos_array[i], unique_res_array[j])]) / N
                    if p > 0:
                        log_prob = round(math.log(20 * p, 2), 2)
                    else:
                        log_prob = float('nan')
                    Z[j, i] = log_prob
                except KeyError:
                    continue
                try:
                    if wt_dict[(unique_pos_array[i], unique_res_array[j])]:
                        xbox.append(i)
                        ybox.append(j)
                except KeyError:
                    continue
    else:
        cmap_color = 'coolwarm'
        for i in xi:
            for j in yi:
                try:
                    Z[j, i] = mut_dict[(unique_pos_array[i], unique_res_array[j])]
                    if wt_dict[(unique_pos_array[i], unique_res_array[j])]:
                        xbox.append(i)
                        ybox.append(j)
                except KeyError:
                    continue
    if resfile:
        positions = []
        with open(resfile, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line:
                    first_split = line.split()[0]
                    if first_split.isdigit():
                        positions.append(int(first_split))
        x_label = []
        for p, position in enumerate(positions):
            if p + 1 in unique_pos_array:
                x_label.append(position)
    else:
        x_label = list(unique_pos_array)
    y_label = list(unique_res_array)

    fig, ax = plt.subplots()
    heat_map = ax.pcolormesh(Z, cmap=cmap_color)
    ax.set_xticks(xi + 0.5, minor=False)
    ax.set_yticks(yi + 0.5, minor=False)
    ax.set_xticklabels(x_label, minor=False)
    ax.set_yticklabels(y_label, minor=False)
    ax.tick_params(direction='out', labelsize=10, top=False, right=False)

    for k in range(0, len(xbox)):
        ax.add_patch(Rectangle((xbox[k], ybox[k]), 1, 1, fill=False, edgecolor=almost_black, lw=2))
    plt.colorbar(heat_map)
    plt.savefig(name, format='png', dpi=1000)


def gen_stats(mut_per_struct_list):
    """calculate avg and var for output in console"""
    avg_mut = np.average(mut_per_struct_list)
    var_mut = np.var(mut_per_struct_list)
    print "Average number of mutations per structure is {} with variance {}".format(avg_mut, var_mut)


def main(wt_fasta, query_fasta_list, heat_map_bool, pwm_bool, resfile, name):
    """calls child functions"""
    mut_dict, mut_per_struct_list, wt_seq = compare_seq(wt_fasta, query_fasta_list)
    gen_stats(mut_per_struct_list)
    if heat_map_bool:
        mk_heat_map(mut_dict, wt_seq, query_fasta_list, pwm_bool, resfile, name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculates number of mutations per seq and generates visualisation")
    parser.add_argument("-m", "--heat_map", action="store_true",
                        help="invokes heat map visualisation")
    parser.add_argument("-pwm", action='store_true',
                        help='display heatmap as position weight matrix. default displays raw counts.')
    parser.add_argument("-n", "--name", default='heat_map.png',
                        help='name of output png (default is heat_map.png')
    parser.add_argument("-r", "--resfile",
                        help='if resfile is provided, labels on heatmap will be positions from resfile')
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-w", "--wt_fasta", required=True,
                           help="fasta file for wild type protein seq")
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more inquiry fasta")

    args = parser.parse_args()
    main(args.wt_fasta, args.fasta, args.heat_map, args.pwm, args.resfile, args.name)
