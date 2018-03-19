#!/usr/bin/env python
import argparse
import gzip
import itertools
import sys
from CodonUsage import sorted_codon_table


def parse_fasta_to_list(fasta_file):
    allowed_aa = {'A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y'}
    sequence_list = []
    with gzip.open(fasta_file, 'r') as f:
        for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
            sequence = sequence.rstrip()
            if set(sequence) <= allowed_aa:
                sequence_list.append(sequence)
    return sequence_list


def protein_to_dna(protein_sequence):
    dna_codons = []
    for amino_acid in protein_sequence:
        codon = sorted_codon_table[amino_acid][0]
        dna_codons.append(codon)
    dna_sequence = ''.join(dna_codons)
    return dna_sequence


def sequence_distance(string1, string2):
    if len(string1) != len(string2):
        raise ValueError("Sequences not of unequal length")

    difference = 0
    for s1, s2 in zip(string1, string2):
        if s1 != s2:
            difference += 1
        if difference > 1:
            break
    if difference == 1:
        return 1
    elif difference == 0:
        return 0
    else:
        return float('NaN')


def protein_single_nt_change(fasta_file):
    seq_list = parse_fasta_to_list(fasta_file)
    print len(seq_list)
    sys.stdout.write('source\ttarget\tinteraction\tdata\n')
    for seq1, seq2 in itertools.combinations(seq_list, 2):
        dna1 = protein_to_dna(seq1)
        dna2 = protein_to_dna(seq2)
        length1 = len(dna1)
        length2 = len(dna2)
        if length1 == length2:
            diff = sequence_distance(dna1, dna2)
            if diff == 1:
                sys.stdout.write('{0}\t{1}\t{2}\t{3}\n'.format(seq1, seq2, 'pp', 'missense'))
        else:
            if length2 < length1:
                length1, length2 = length2, length1
                dna1, dna2 = dna2, dna1
                seq1, seq2 = seq2, seq1
            partial_diff = sequence_distance(dna1, dna2[:length1])
            if partial_diff == 0:
                extra_codon = dna2[length1:length1 + 3]
                convert_to_stop = False
                for stop_codon in sorted_codon_table['*']:
                    difference = sum(el1 != el2 for el1, el2 in zip(extra_codon, stop_codon))
                    if difference == 1:
                        convert_to_stop = True
                if convert_to_stop:
                    sys.stdout.write('{0}\t{1}\t{2}\t{3}\n'.format(seq1, seq2, 'pp', 'nonsense'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to determine protein sequences present in a fasta fasta
    file that can be converted into another sequence also present in the fasta by changing one nucleotide.
    output formatted for input into cytoscape""")
    required = parser.add_argument_group('required')
    required = parser.add_argument('-f', '--fasta', required=True,
                                   help='fasta file containing protein sequences to be compared')
    args = parser.parse_args()
    protein_single_nt_change(args.fasta)
