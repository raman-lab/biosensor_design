#!/usr/bin/env python

import argparse
import itertools
import sys


def thin_fasta(fasta, l, a, p, c, b, positions):
    if not positions:
        for fasta_file in fasta:
            with open(fasta_file, 'r') as f:
                for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
                    if (
                            (not l or l in tag) and
                            (not a or a in tag) and
                            (not p or p in tag) and
                            (not c or c in tag) and
                            (not b or b in tag)
                    ):
                        sys.stdout.write(tag)
                        sys.stdout.write(sequence)
    else:
        for fasta_file in fasta:
            with open(fasta_file, 'r') as f:
                for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
                    seq_list = list(sequence)
                    for i, aa in enumerate(seq_list):
                        if i not in positions:
                            seq_list[i] = '-'
                    sequence = ''.join(seq_list)
                    sys.stdout.write(tag)
                    sys.stdout.write('{0}\n'.format(sequence))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to thin a fasta file be specifiying sequences qualifiers in the sequence
        id. for example, A1, B0, or cbio. qualifiers are specified through options, which all accept a
        one input (ie -A A0. not -A A0 A1). thinned sequences must pass all qualifiers and are written to stdout"""
    )
    parser.add_argument('-b', help='specify seqs from certain block to keep')
    parser.add_argument('-c', help='specify seqs from certain conformation to keep')
    parser.add_argument('-p', help='specify seqs from certain pairing to keep')
    parser.add_argument('-a', help='specify seqs from certain alignment to keep')
    parser.add_argument('-l', help='specify seqs from certain ligand to keep')
    parser.add_argument('-positions', nargs='*', type=int,
                        help='takes space delimited integers. output will only be those postions. '
                             'all other postions will be dashes')
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more fasta files")

    args = parser.parse_args()
    thin_fasta(args.fasta, args.l, args.a, args.p, args.c, args.b, args.positions)
