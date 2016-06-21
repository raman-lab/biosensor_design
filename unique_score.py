#!/usr/bin/env python

"""take a unique fasta files and a composite score file
return a unique score file with only entries specified in fastas"""

import argparse


def get_tags(fastas):
    """extract pdb identifier from fasta files"""
    tags = []
    for fasta in fastas:
        with open(fasta, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    if '/' in line:
                        pre_tag = line.split('/')[-1]
                    else:
                        pre_tag = line.split('>')[-1]
                    tag = pre_tag.split('.')[0]
                    tags.append(tag)
    return tags


def unique_sc(score, fastas):
    """use tags to populate a unique score file from the composite score file"""
    tags = get_tags(fastas)
    with open(score, 'r') as f:
        out_file = score.split('.')[0] + '_unique.' + score.split('.')[1]
        o = open(out_file, 'w')
        top = f.readline()
        o.write(top)
        header = top.split()
        d_index = header.index('description')
        for line in f:
            if (not line) or (line.startswith("#")) or (line[0].isalpha()):
                continue
            line_list = line.split()
            if line_list[d_index] in tags:
                o.write(line)
        o.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="parse a composite score file and return entries that match input fastas")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="pdb for wild type protein")
    requiredO.add_argument("-s", "--score", required=True, help='composite score file')

    args = parser.parse_args()
    unique_sc(args.score, args.fasta)