#!/usr/bin/env python

import argparse
from itertools import izip_longest
import sys
from Bio import Seq


def check_fragments(oligo_file, design_fasta):
    design_aa_list = []
    with open(design_fasta, 'r') as f:
        for pdb, seq in izip_longest(f, f, fillvalue=None):
            if '4AC0' and 'B0' in pdb:
                block = seq[77:117]
            elif '4AC0' and 'B1' in pdb:
                block = seq[99:138]
            elif '2uxo' and 'B0' in pdb:
                block = seq[62:100]
            elif '2uxo' and 'B1' in pdb:
                block = seq[136:176]
            else:
                raise Exception('Unrecognized design name')
            design_aa_list.append(block)

    fragment_list = []
    with open(oligo_file, 'r') as o:
        for pdb, seq in izip_longest(o, o, fillvalue=None):
            if '4AC0' and 'B0' in pdb:
                seq_lower = seq.lower()
                seq_no_5p = seq_lower.split('gtgacccgtccctgggtctcaagat')[1]
                fragment = seq_no_5p.split('gccttgagaccgggcagaggtcgac')[0]
            elif '4AC0' and 'B1' in pdb:
                seq_lower = seq.lower()
                seq_no_5p = seq_lower.split('tgcccgctgtcttcaggtctcaagta')[1]
                fragment = seq_no_5p.split('catttgagacctgtagcccggcagtg')[0]
            elif '2uxo' and 'B0' in pdb:
                seq_lower = seq.lower()
                seq_no_5p = seq_lower.split('cgatcgtgcccacctggtctccactg')[1]
                fragment = seq_no_5p.split('gttctgagaccagttggagcccgcac')[0]
            elif '2uxo' and 'B1' in pdb:
                seq_lower = seq.lower()
                seq_no_5p = seq_lower.split('ctggtgcgtcgtctggtctctggat')[1]
                fragment = seq_no_5p.split('cgttggagaccggcgaacacttccc')[0]
            else:
                raise Exception('Unrecognized oligo name')
            fragment_list.append(fragment)

    missing_list = []
    for item in fragment_list:
        aa_fragment = Seq.translate(item)
        if aa_fragment in design_aa_list:
            design_aa_list.remove(aa_fragment)
        else:
            missing_list.append(aa_fragment)
    if missing_list:
        sys.stderr.write('Error: The following oligo sequences do not match a design amino acid sequence\n')
        for miss in missing_list:
            sys.stderr.write('{0}\n'.format(miss))
    if design_aa_list:
        sys.stderr.write('Error: The following design sequences do not match an oligo sequence\n')
        for design in design_aa_list:
            sys.stderr.write('{0}\n'.format(design))
    sys.stdout.write('done\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to check Kyle's mutateproteinintodnapy2.py
        script checks equivalence of translated oligo fragment and designed amino acid block
        an error will be printed to std.err if:
            the oligo does not have the correct primers, BsaI cut sites, or 5 nucleotide buffers
            oligo sequences and design sequences are not one to one
            the two sequences are not equivalent"""
    )
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-o", "--oligo", required=True,
                           help="Oligo output file")
    requiredO.add_argument("-d", "--design", required=True,
                           help="Amino acid fasta file")

    args = parser.parse_args()
    check_fragments(args.oligo, args.design)
