#!/usr/bin/env python

# Automates biosensor design post processing from rosetta output

import argparse
import gen_enzdes_cutoffs
import postProcessPlot


def main(inputScFile, descriptorFile, totalOutputSeqNum, plot):
    print 'Making Cutoffs'
    gen_enzdes_cutoffs.main(inputScFile, descriptorFile, 'cutoffs.cut')
    # 
    # execute perl script to make loosely filtered score files
    # revert to native script
    # unique seq check
    # execute perl script to tightly filter score files

    if plot:
        print 'Plotting Data'
        postProcessPlot.main('total_score', 'SR_1_total_score', 'filteredComposite.sc', 'inputScFile')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Automate biosensor design post processing from rosetta \
            output.
            From pdbs and scorefiles, final output will be DNA sequences with \
            flanking primer sequences for synthesis.
            Intermediate output includes: scatter plots of designed ligand \
            alignments (filtered and unfilered), cut file, unique filtered pdbs and \
            fastas, amino acid frequencies, coovariation graph, and weblogo""")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="turn on plotting")
    parser.add_argument("-t", "--total", type=int,
                        help="number of output protein and DNA sequences (default: 10000)",
                        default=10000)

    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-c", "--cutoff", required=True, help="descriptor specification for cutoff")
    requiredO.add_argument("-s", "--scoreFile", nargs='*', required=True,
                           help="composite score file. To generate run: find . -name '*.sc' | xargs cat > compScore.sc")
    args = parser.parse_args()

    main(args.scoreFile, args.cutoff, args.total, args.plot)
