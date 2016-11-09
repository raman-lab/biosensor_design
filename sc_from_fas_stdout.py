#!/usr/bin/env python

"""take a unique fasta files and a composite score file
prints score file lines to stdout for redirection"""

import argparse
import sys


def get_tags(fastas):
    """extract pdb identifier from fasta files"""
    tags = set()
    for fasta in fastas:
        with open(fasta, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if '/' in line:
                        pre_tag = line.split('/')[-1]
                    else:
                        pre_tag = line.split('>')[-1]
                    tag = pre_tag.split('.')[0]
                    tags.add(tag)
    return tags


def unique_sc(score, fasta_list):
    """use tags to print score file lines from the composite score file"""
    tags = get_tags(fasta_list)
    with open(score, 'r') as f:
        top = f.readline()
        sys.stdout.write(top)
        header = top.split()
        d_index = header.index('description')
        for line in f:
            if (not line) or (line.startswith("#")) or (line[0].isalpha()):
                continue
            line_list = line.split()
            desc = line_list[d_index].split('.')[0]
            if desc in tags:
                tags.discard(desc)
                sys.stdout.write(line)
        if len(tags) > 0:
            sys.stderr.write('Error: Not all sequences have a matching score term\n')
            sys.stderr.write('Add the score terms for structures listed below to the specified score file:\n')
            for tag in tags:
                sys.stderr.write('{0}\n'.format(tag))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="prints score file lines to stdout for redirection "
                                                 "takes a unique fasta files and a composite score file")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more input fasta, which specify desired structures in output score file")
    requiredO.add_argument("-s", "--score", required=True, help='composite score file')

    args = parser.parse_args()
    unique_sc(args.score, args.fasta)
