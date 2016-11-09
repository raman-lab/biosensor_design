#!/usr/bin/env python
import argparse

import itertools

import sys


def get_indices(resfile):
    mutable_residues_indices = []
    with open(resfile, 'r') as f:
        for line in f:
            if line[0].isdigit():
                residue_index = int(line.split()[0]) - 1
                mutable_residues_indices.append(residue_index)
    return mutable_residues_indices


def get_mutable_set(fasta_files, resfiles, block, mut_bool):
    for fasta in fasta_files:
        for resfile in resfiles:
            mutable_residue_indices = get_indices(resfile)
            with open(fasta, 'r') as f:
                for id, sequence in itertools.izip_longest(f, f, fillvalue=None):
                    mutable_residues = []
                    if block and block not in id:
                        continue
                    if mut_bool:
                        for index in mutable_residue_indices:
                            mutable_residues.append(sequence[index])
                        mutable_set = ''.join(mutable_residues)
                    else:
                        mutable_set = sequence[mutable_residue_indices[0]:mutable_residue_indices[-1]+1]
                    sys.stdout.write('{0}'.format(id))
                    sys.stdout.write('{0}\n'.format(mutable_set))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to reduce a full rosetta design sequence to just the residues that were considered
        for design. script was developed specifically to generate sequences sets for use in ml algorithms.
        writes to stdout"""
    )
    parser.add_argument('-b', '--block', help='restricts output to B0 or B1. only B0 and B1 can be used as inputs')
    parser.add_argument('-mut', '--mutable', action="store_true",
                        help='restrict output to only the residues that were considered for mutation'
                             ' in rosetta design runs')
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more fasta files")
    requiredO.add_argument("-r", "--resfile", nargs='*', required=True,
                           help='resfile used in rosetta during design process')

    args = parser.parse_args()
    get_mutable_set(args.fasta, args.resfile, args.block, args.mutable)
