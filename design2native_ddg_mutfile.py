#!/usr/bin/env python 
import argparse
import collections
import itertools


def get_wt_sequence(wt_fasta):
    with open(wt_fasta, 'r') as f:
        f.readline()
        wt_sequence = f.readline()
    return wt_sequence


def compare_sequences(wt_sequence, sequence):
    mutation_dict = collections.OrderedDict()
    for index, residue in enumerate(sequence):
        if residue is not wt_sequence[index]:
            mutation_dict[index] = (wt_sequence[index], residue)
    return mutation_dict


def write_mutfile(index_residue_dict, identifier):
    mutation_count = len(index_residue_dict.keys())
    lines_to_write = ['total 1\n', '{0}\n'.format(mutation_count)]
    for index, residue_tuple in index_residue_dict.items():
        mutation_line = '{0} {1} {2}\n'.format(residue_tuple[0], index + 1, residue_tuple[1])
        lines_to_write.append(mutation_line)
    with open('./{0}.mutfile'.format(identifier.split('/')[-1].rstrip(), 'w')) as f:
        f.writelines(lines_to_write)


def make_design2native_ddg_mutfile(wt_fasta, fastas):
    wt_sequence = get_wt_sequence(wt_fasta)
    for fasta in fastas:
        with open(fasta, 'r') as f:
            for identifier, sequence in itertools.izip_longest(f, f, fillvalue=None):
                index_residue_dict = compare_sequences(wt_sequence, sequence)
                write_mutfile(index_residue_dict, identifier)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compare a fasta of designed protein sequences to the native sequence."
                                                 "For each sequence, each mutated residue will be added to a mutfile "
                                                 "for use in Rosetta's ddg_monomer application.")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-w", "--wt_fasta", required=True,
                           help="fasta file with wild type protein seq. Note: only first sequence in file will be used")
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more fasta files with designed sequences")

    args = parser.parse_args()
    make_design2native_ddg_mutfile(args.wt_fasta, args.fasta)
