#!/usr/bin/env python
import argparse
import sys


def renumber_pdb(pdb_file, last_residue_number):
    seen_atom_1 = 0
    with open(pdb_file, 'r') as f:
        for line in f:
            split_line = line.split()
            if line == 'END\n':
                break
            atom_number = int(split_line[1])
            if atom_number == 1:
                seen_atom_1 += 1
            if seen_atom_1 == 2:
                number = int(split_line[5])
                line = line.replace('A{:4d}'.format(number), 'A{:4d}'.format(number + last_residue_number))
            sys.stdout.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""quick script to renumber residues in pdb file. input is a pdb
    file (homodimer) with both monomers on the same chain with the same numbering and the last residue of number of
    the monomer. output is written to stdout.""")
    required = parser.add_argument_group("required")
    parser.add_argument('-p', '--pdb', required=True, help="name of pdb file to be renumbered")
    parser.add_argument('-l', '--last_residue_number', required=True, type=int,
                        help="last residue number of a single monomer of the protein")
    args = parser.parse_args()
    renumber_pdb(args.pdb, args.last_residue_number)
