#!/usr/bin/env python

import argparse
import sys
from Bio import AlignIO
from Bio.Align import AlignInfo


def consensus_main(msa_files):
    for msa in msa_files:
        alignment = AlignIO.read(msa, 'clustal')
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus = summary_align.dumb_consensus(threshold=0.5)
        str_consensus = str(consensus)
        dash_consensus = str_consensus.replace('X', '-')
        sys.stdout.write('{0}\n'.format(dash_consensus))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""first attempt at making a consensus seq in python with msa as input.
        currently only accepts clustal files."""
    )
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-m", "--msa", nargs='*', required=True,
                           help="one or more clustal msa files")

    args = parser.parse_args()
    consensus_main(args.msa)