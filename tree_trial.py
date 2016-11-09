#!/usr/bin/env python

import argparse
from Bio.Phylo import TreeConstruction
from Bio import AlignIO
from Bio import Phylo


def tree_main(msa_files):
    for msa in msa_files:
        alignment = AlignIO.read(msa, 'clustal')
        calculator = TreeConstruction.DistanceCalculator('identity')
        dist_matrix = calculator.get_distance(alignment)
        constructor = TreeConstruction.DistanceTreeConstructor(calculator, 'nj')
        tree = constructor.build_tree(alignment)
        tree.rooted = True
        Phylo.draw(tree)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""first attempt at making a phylogenetic tree in python with msa as input.
        currently only accepts clustal files."""
    )
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-m", "--msa", nargs='*', required=True,
                           help="one or more clustal msa files")

    args = parser.parse_args()
    tree_main(args.msa)
