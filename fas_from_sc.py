#!/usr/bin/env python
import argparse
import itertools

import sys


def get_descriptions_from_score_file(score_file):
    descriptions = set()
    with open(score_file, 'r') as f:
        header_list = f.readline().split()
        d_index = header_list.index('description')
        for line in f:
            if (not line) or (line.startswith("#")) or (line[0].isalpha()):
                continue
            line_list = line.split()
            desc = line_list[d_index].split('.')[0]
            descriptions.add(desc)
    return descriptions


def fas_from_sc(fastas, score_files):
    for score_file in score_files:
        descriptions = get_descriptions_from_score_file(score_file)
        for fasta in fastas:
            with open(fasta, 'r') as f:
                for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
                    if '/' in tag:
                        pre_tag = tag.split('/')[-1]
                    else:
                        pre_tag = tag.split('>')[-1]
                    new_tag = pre_tag.split('.')[0]

                    if new_tag in descriptions:
                        descriptions.discard(new_tag)
                        sys.stdout.write('>{0}\n'.format(new_tag))
                        sys.stdout.write(sequence)

                if len(descriptions) > 0:
                    sys.stderr.write('Error: Not all descriptions in score file have a matching fasta sequence\n')
                    sys.stderr.write('The following descriptions do not have a match:\n')
                    for desc in descriptions:
                        sys.stderr.write('{0}\n'.format(desc))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""script that returns fasta sequences that correspond to entries in a Rosetta score file."""
    )
    required = parser.add_argument_group('required arguments')
    required.add_argument('-f', '--fastas', nargs='*', required=True,
                          help='composite fasta(s) containing protein sequences named in score files')
    required.add_argument('-s', '--score_file', nargs='*', required=True,
                          help='Rosetta score files that contain name of fasta sequences to be pulled '
                               'from input fasta files')
    args = parser.parse_args()
    fas_from_sc(args.fastas, args.score_file)