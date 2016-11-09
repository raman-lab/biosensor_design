#!/usr/bin/env python
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Script to demultiplex next generation sequencing data based on input barcodes
        and input constant regions"""
    )
    parser.add_argument("-c", "--constant", type=int, default=100,
                        help="min number of seqs for a ligand per A<int>_P<int> combo")
    parser.add_argument("-max", "--max_seq", type=int, default=50000,
                        help="total size of curated set for the scaffold (not by ligand or combo)"
                             "before minimum is added")
    parser.add_argument("-b", "--blind", action="store_true",
                        help="blindly add minimum to each A<int>_P<int> combo and do not consider how many designs"
                             "are already present in the combo")
    requiredO = parser.add_argument_group('required arguments')