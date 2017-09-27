#!/usr/bin/env python
import argparse
import collections
import itertools
import numpy as np
import pandas as pd
import re
import sys


class Bin(object):
    def __init__(self, barcode_sequence_str, median_fluor_float, fraction_float):
        self._sequence = None
        self._median_fluor = None
        self._fraction = None
        self.sequence = barcode_sequence_str
        self.median_fluor = float(median_fluor_float)
        self.fraction = float(fraction_float)
        # self._protein = None
        # self._state = None
        # self._replicate = None

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence_str):
        self._sequence = str(sequence_str).upper()

    @property
    def median_fluor(self):
        return self._median_fluor

    @median_fluor.setter
    def median_fluor(self, median_fluor_float):
        self._median_fluor = float(median_fluor_float)

    @property
    def fraction(self):
        return self._fraction

    @fraction.setter
    def fraction(self, fraction_float):
        if fraction_float > 1:
            self._fraction = float(fraction_float) / 100
        else:
            self._fraction = float(fraction_float)

    # @property
    # def protein(self):
    #     return self._protein
    #
    # @protein.setter
    # def protein(self, protein_str):
    #     self._protein = str(protein_str)
    #
    # @property
    # def state(self):
    #     return self._state
    #
    # @state.setter
    # def state(self, state_str):
    #     self._state = str(state_str)
    #
    # @property
    # def replicate(self):
    #     return self._replicate
    #
    # @replicate.setter
    # def replicate(self, replicate_int):
    #     self._replicate = int(replicate_int)


def parse_barcode_file(barcode_files):
    label_dict = collections.OrderedDict(
        [('barcode', None), ('protein', None), ('state', None), ('fluor', None), ('percent', None), ('licate', None)])
    protein_replicate_state_dict = collections.defaultdict(lambda: collections.defaultdict(dict))
    for barcode_file in barcode_files:
        with open(barcode_file, 'r') as f:
            first_line_list = f.next().rstrip().lower().split('\t')
            labels = label_dict.keys()
            for index, item in enumerate(first_line_list):
                for label in labels:
                    if label in item:
                        label_dict[label] = index
                        labels.remove(label)
                        break
            if None in label_dict.values():
                raise Exception('One of the essential labels in the barcode input file was not found')
            for line in f:
                if not line.isspace() and not line.startswith('#'):
                    line_list = line.rstrip().split('\t')
                    barcode_seq = line_list[label_dict['barcode']].upper()
                    protein = line_list[label_dict['protein']]
                    state = line_list[label_dict['state']].lower()
                    median_fluor = line_list[label_dict['fluor']]
                    percent = line_list[label_dict['percent']]
                    replicate = int(line_list[label_dict['licate']])

                    if 'rep' in state:
                        state = 'repressed'
                    elif 'ind' in state:
                        state = 'induced'
                    elif 'free' in state:
                        state = 'free'
                    else:
                        raise Exception('States can only be repressed, induced, or free')
                    # for item in line_list:
                    #
                    # barcode_seq, protein, state, median_fluor, percent, replicate = line.rstrip().split()
                    # replicate = int(replicate)

                    bin_obj = Bin(barcode_seq, median_fluor, percent)

                    if state not in protein_replicate_state_dict[protein][replicate].keys():
                        protein_replicate_state_dict[protein][replicate][state] = []
                    protein_replicate_state_dict[protein][replicate][state].append(bin_obj)

    return protein_replicate_state_dict


def quality_check_fastq(quality_str, max_expected_incorrect):
    expected_incorrect_num = 0
    for ascii_char in quality_str.rstrip():
        phred_score = ord(ascii_char) - 33
        prob_incorrect = pow(10, - float(phred_score) / 10)
        # maybe add an abs prob cutoff (ie every base must be < 0.01). maybe only do this for the library and not
        # whole read
        expected_incorrect_num += prob_incorrect
    if expected_incorrect_num > max_expected_incorrect:
        return False
    else:
        return True


def check_barcode(full_sequence, barcode_dict):
    for barcode in barcode_dict.keys():
        if full_sequence.startswith(barcode) or full_sequence.endswith(barcode):
            return barcode
    return None


def parse_from_library_architecture(library_architecture):
    with open(library_architecture, 'r') as f:
        library_read = f.readline().rstrip().upper()
    library_start = library_read.find('N')
    library_end = library_read.rfind('N')
    barcode_start = library_read.find('X')
    barcode_end = library_read.rfind('X')

    lib_start_str = library_read[library_start - 6:library_start]
    lib_end_str = library_read[library_end + 1: library_end + 7]
    barcode_start_str = library_read[barcode_start - 6: barcode_start]
    barcode_end_str = library_read[barcode_end + 1: barcode_end + 7]

    return (lib_start_str, lib_end_str), (barcode_start_str, barcode_end_str)


def calculate_genotype_fluorescence(fastq_files, barcode_files, library_architecture,
                                    max_expected_incorrect, minimum_count):
    protein_replicate_state_dict = parse_barcode_file(barcode_files)
    const_library_tup, const_barcode_tup = parse_from_library_architecture(library_architecture)
    sequence_dict = collections.defaultdict(dict)

    # variables for my own use - will be removed in future
    total_sequence_count = 0
    sequences_passing_qc = 0
    sequences_with_matching_const_region = 0

    for fastq_file in fastq_files:
        with open(fastq_file, 'r') as f:
            for identifier, sequence, spacer, quality_str in itertools.izip_longest(*[f] * 4, fillvalue=None):
                total_sequence_count += 1
                if quality_check_fastq(quality_str, max_expected_incorrect):
                    library_seq = re.search('{0}(.*){1}'.format(const_library_tup[0], const_library_tup[1]), sequence)
                    barcode_seq = re.search('{0}(.*){1}'.format(const_barcode_tup[0], const_barcode_tup[1]), sequence)
                    sequences_passing_qc += 1
                    if library_seq and barcode_seq:
                        library_seq = library_seq.group(1)
                        barcode_seq = barcode_seq.group(1)
                        if barcode_seq not in sequence_dict[library_seq].keys():
                            sequence_dict[library_seq] = collections.Counter()
                        sequence_dict[library_seq][barcode_seq] += 1
                        sequences_with_matching_const_region += 1

    # using first normalization scheme
    # count_matrix = pd.DataFrame.from_dict(sequence_dict, orient='index') # c_ij
    # barcode_counts = count_matrix.sum(axis='index')
    # barcode_fractions = barcode_counts / sum(barcode_counts) # f_j
    # fraction_matrix = barcode_fractions * count_matrix.div(count_matrix.sum(axis='index'), axis='columns') # b_ij
    # normalized_fraction_matrix = fraction_matrix.div(fraction_matrix.sum(axis='columns'), axis='index') # a_ij
    #
    # barcode_fluor_dict = dict((k, v) for k, v in barcode_dict.items() if k in normalized_fraction_matrix.columns.values)
    # median_fluor = pd.Series.from_array(barcode_fluor_dict.values(), index=barcode_fluor_dict.keys()) # m_j
    #
    # log_fluor_matrix = normalized_fraction_matrix * np.log10(median_fluor)
    # seq_vec = log_fluor_matrix.sum(axis='columns')
    # mean_fluor = np.power(10, seq_vec)
    #
    # sys.stdout.write('{}\n'.format(mean_fluor))
    # sys.stdout.write('total seqs: {}\n'.format(total_sequence_count))
    # sys.stdout.write('missed: {}\n'.format(sequences_without_barcode_count))

    count_matrix = pd.DataFrame.from_dict(sequence_dict, orient='index', dtype=np.uint32)
    count_matrix[count_matrix < minimum_count] = 0
    # count_matrix[count_matrix >= minimum_count] = 1

    # barcode_fluor_dict = {}
    # barcode_fraction_dict = {}
    # barcode_rename_dict = {}
    # for barcode in barcode_list:
    #     barcode_fluor_dict[barcode.sequence] = barcode.median_fluor
    #     barcode_fraction_dict[barcode.sequence] = barcode.bin_fraction
    #     barcode_rename_dict[barcode.sequence] = '{0}_{1}'.format(barcode.state, barcode.replicate)

    for protein in protein_replicate_state_dict.keys():
        replicate_matrices = []
        replicate_keys = []
        for replicate in protein_replicate_state_dict[protein].keys():
            replicate_keys.append('replicate {0}'.format(replicate))
            state_matrices = []
            rename_dict = {}
            for state in protein_replicate_state_dict[protein][replicate].keys():
                state_fractions = []
                state_median_fluors = []
                state_seqs = []

                for bin_obj in protein_replicate_state_dict[protein][replicate][state]:
                    state_fractions.append(bin_obj.fraction)
                    state_median_fluors.append(bin_obj.median_fluor)
                    state_seqs.append(bin_obj.sequence)
                    rename_dict[bin_obj.sequence] = state
                state_median_fluor_vec = pd.Series(state_median_fluors, index=state_seqs)
                state_fraction_vec = pd.Series(state_fractions, index=state_seqs)
                state_fraction_vec = state_fraction_vec.div(state_fraction_vec.sum())
                state_matrix = count_matrix.loc[:, state_seqs]
                state_matrix.dropna(axis='index', how='all', inplace=True)
                dummy = state_matrix.div(state_matrix.sum(axis='index'), axis='columns')
                fraction_matrix = state_fraction_vec * state_matrix.div(state_matrix.sum(axis='index'), axis='columns')
                norm_fraction_matrix = fraction_matrix.div(fraction_matrix.sum(axis='columns'), axis='index')
                log_fluor_matrix = norm_fraction_matrix * np.log10(state_median_fluor_vec)
                seq_vec = log_fluor_matrix.sum(axis='columns')
                state_mean_fluor = np.power(10, seq_vec)
                state_mean_fluor.rename(state, inplace=True)
                state_matrices.append(state_mean_fluor)

            replicate_matrix = pd.concat(state_matrices, axis='columns')
            replicate_matrix['fold_induction'] = replicate_matrix['induced'] / replicate_matrix['repressed']
            replicate_matrices.append(replicate_matrix)
        protein_matrix = pd.concat(replicate_matrices, axis='columns', keys=replicate_keys)
        sort_tup = zip(replicate_keys, ['fold_induction'] * len(replicate_keys))
        protein_matrix.sort_values(by=sort_tup, ascending=[True] * len(replicate_keys), inplace=True)
        sys.stdout.write(protein)
        sys.stdout.write('{0}\n'.format(protein_matrix.to_string()))
        print count_matrix.info()
        # excel_writer = pd.ExcelWriter('test.xlsx')
        # protein_matrix.to_excel(excel_writer)
        # excel_writer.save()
        ## left off here
        # protein_matrix.sort([replicate_keys[]], axis='columns')



    # fraction_vec = pd.Series.from_array(barcode_fraction_dict.values(), index=barcode_fraction_dict.keys())  # f_j
    # median_fluor_vec = pd.Series.from_array(barcode_fluor_dict.values(), index=barcode_fraction_dict.keys())  # m_j
    #
    # fraction_matrix = fraction_vec * count_matrix.div(count_matrix.sum(axis='index'), axis='columns')  # b_ij
    # normalized_fraction_matrix = fraction_matrix.div(fraction_matrix.sum(axis='columns'), axis='index')  # a_ij
    #
    # fluor_matrix = normalized_fraction_matrix * median_fluor_vec
    # fluor_matrix.rename(columns=barcode_rename_dict, inplace=True)

    # need to insert colms for fold values now


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to calculate the mean fluorescence value for a specific genotype present in sequence data.
        input is fastq(s) of sequences and a tab delimited text file of barcodes, median fluorescence values,
        population percentages, and replicate numbers for each barcode."""
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fastq', nargs='*', required=True, help='text file(s) of sequences')
    required.add_argument('-b', '--barcodes', nargs='*', required=True,
                          help='tab delimited text file with one column for barcodes, protein, states, '
                               'median fluorescence values, population percentages, and replicate number.')
    required.add_argument('-l', '--library_architecture', required=True,
                          help='text file with an example one-line library read used to specify sequence architecture.'
                               'library sequence should be denoted with "N"s, and barcodes should be denoted '
                               'with "X"s. currently this script is takes six residues before and after the library '
                               'and barcode sequences, assumed constant, to extract library and barcode sequences.'
                               'more functionality will be added as needed.')
    parser.add_argument('-e', '--expected_error_max', default=1, type=int,
                        help='integer number of the number of expected wrong bases per read based on phred quality '
                             'score')
    parser.add_argument('-m', '--minimum', default=5, type=int,
                        help='minimum number of reads a sequence must have to be included in fluorescence calculations')
    args = parser.parse_args()
    calculate_genotype_fluorescence(args.fastq, args.barcodes, args.library_architecture,
                                    args.expected_error_max, args.minimum)
