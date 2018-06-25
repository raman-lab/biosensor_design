#!/usr/bin/env python
import argparse
import itertools
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def parse_first_fasta_sequence(fasta_file):
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    return lines[1].rstrip()


def parse_fasta_to_tag_seq_dict(fasta_file):
    tag_sequence_dict = {}
    with open(fasta_file, 'r') as f:
        for tag, sequence in itertools.izip_longest(f, f, fillvalue=None):
            tag = tag.rstrip().split('>', 1)[-1].split(' ')[0]
            sequence = sequence.rstrip()
            tag_sequence_dict[tag] = sequence
    return tag_sequence_dict


def make_rosetta_dna_jobs(fasta, pdb_dna_fasta, pdb, start_pos, dna_chain):
    pdb_dna_seq = parse_first_fasta_sequence(pdb_dna_fasta)
    tag_seq_dict = parse_fasta_to_tag_seq_dict(fasta)
    dna_abbreviation_dict = {'A': 'ADE', 'C': 'CYT', 'G': 'GUA', 'T': 'THY'}

    for tag, seq in tag_seq_dict.iteritems():
        shell_commands = []
        if len(seq) < len(pdb_dna_seq):
            raise Exception('each sequence must be at least as long as the pdb dna sequence')
        len_diff = len(pdb_dna_seq) - len(seq)
        rev_seq = str(Seq(seq, alphabet=generic_dna).reverse_complement())
        for i in range(0, len_diff):
            for s in [seq, rev_seq]:
                if s == seq:
                    strand = 'fwd'
                else:
                    strand = 'rev'
                dna_def_list = []
                for index, base_tup in enumerate(zip(pdb_dna_seq[i:], s)):
                    if base_tup[0] != base_tup[1]:
                        position = start_pos + index
                        dna_def_str = '{0}.{1}.{2}'.format(dna_chain, position, dna_abbreviation_dict[base_tup[1]])
                        dna_def_list.append(dna_def_str)

                log_file = '{0}_tag{1}_{2}_frame{3}.log'.format(pdb.split('/')[-1].split('.')[0], tag, strand, i)
                shell_command = '/scratch/nwhoppe/Rosetta/main/source/bin/rosettaDNA.linuxgccrelease ' \
                                '-database /scratch/nwhoppe/Rosetta/main/database' \
                                '@dna_affinity.flags -in::file::s {0} -out::suffix _tag{1}_{2}_frame{3}' \
                                '-parser::protocol dna_affinity.xml -parser::script_vars' \
                                'dna_seq={4} > {5}\n'.format(pdb, tag, strand, i, ','.join(dna_def_list), log_file)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='script to make condor job files for rosetta dna affinity calculation. '
                    'script is specific to pdb 3lsr and paths on biochemistry computer cluster')
    parser.add_argument('-p', '--pdb', default='3lsr_native.pdb', help='name of pdb structure used as input to rosetta')
    parser.add_argument('-s', '--start_pos', type=int, default=5, help='pdb number of first dna base')
    parser.add_argument('-c', '--dna_chain', default='C', help='chain of top strand of dna in pdb')
    required = parser.add_argument_group('required')
    required.add_argument('-f', '--fasta', required=True,
                          help='fasta file containing sequences to be modeled')
    required.add_argument('-pf', '--pdb_dna_fasta', required=True,
                          help='fasta file with the dna sequence present in the pdb file')
    args = parser.parse_args()
    make_rosetta_dna_jobs(args.fasta, args.pdb_dna_fasta, args.pdb, args.start_pos, args.dna_chain)