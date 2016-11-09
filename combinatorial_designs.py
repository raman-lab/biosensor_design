#!/usr/bin/env python
### UNTESTED
import argparse
import itertools
import sys


def get_wt_seqs(wt_fasta_list):
    wt_dict = {}
    for wt_seq in wt_fasta_list:
        with open(wt_seq, 'r') as f:
            for identifier, sequence in itertools.izip_longest(f, f, fillvalue=None):
                pdb = identifier.split('>')[-1]
                wt_dict[pdb] = sequence
    return wt_dict


def get_variants(identifier, sequence, wt_seq):
    mutation_position_dict = {}
    for index, residue in enumerate(sequence):
        if residue is not wt_seq[index]:
            mutation_position_dict[residue] = index + 1

    mutated_residues = ''.join(mutation_position_dict.keys())
    combinations = []
    for p in range(1, len(mutated_residues) + 1):
        n_combo = [list(x) for x in itertools.combinations(mutated_residues, p)]
        combinations.append(n_combo)

    combinatorial_sequence_dict = {}
    combinatorial_sequence = sequence
    for combo in combinations:
        new_identifier = '{0}_{1}'.format(identifier, len(combo))
        for mutation in combo:
            index = mutation_position_dict[mutation] - 1
            combinatorial_sequence = combinatorial_sequence[index].replace(wt_seq[index], mutation)
            mutation_string = '{0}{1}{2}'.format(wt_seq[index], index + 1, mutation)
            new_identifier = '{0}_{1}'.format(new_identifier, mutation_string)
        combinatorial_sequence_dict[new_identifier] = combinatorial_sequence

    return combinatorial_sequence_dict


def make_combinatorial_variants(fasta_list, wt_fasta_list):
    wt_seqs = get_wt_seqs(wt_fasta_list)
    combinatorial_variants = []
    for pdb, wt_seq in wt_seqs.items():
        for fasta in fasta_list:
            with open(fasta, 'r') as f:
                for identifier, sequence in itertools.izip_longest(f, f, fillvalue=None):
                    if pdb in identifier.split('_'):
                        variant_dict = get_variants(identifier, sequence, wt_seq)
                        for updated_identifier, variant in variant_dict.items():
                            if variant not in combinatorial_variants:
                                combinatorial_variants.append(variant)
                                sys.stdout.write('{0}\n'.format(updated_identifier))
                                sys.stdout.write('{0}\n'.format(variant))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to generate a combinatorial set of mutated protein sequences given input protein
        sequences with at least one mutation in the form of a fasta sequences.
        script checks internally for uniqueness. should be used in conjunction with uniquify_fas.py to check for
        uniqueness externally"""
    )
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more fasta files with mutated protein sequences")
    requiredO.add_argument("-w", "--wt_fasta", nargs='*', required=True,
                           help="fasta file(s) for wild type protein seq. wt sequence should be labeled with pdb "
                                "id (ie >4AC0")

    args = parser.parse_args()
    make_combinatorial_variants(args.fasta, args.wt_fasta)
