#!/usr/bin/env python
"""given input fasta(s), script outputs non-unique sequences in fasta fromat to stdout"""
import sys
import itertools


def non_unique(filename):
    with open(filename, 'r') as f:
        for identifier, sequence in itertools.izip_longest(f, f, fillvalue=None):
            if sequence in seen:
                sys.stdout.writelines([identifier, sequence])
            else:
                seen.add(sequence)


if __name__ == '__main__':
    seen = set()
    for fasta in sys.argv[1:]:
        non_unique(fasta)
