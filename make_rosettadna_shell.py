#!/usr/bin/env python
import argparse
import sys


def parse_paired_alignment(fasta_alignment_pair):
    with open(fasta_alignment_pair, 'r') as f:
        lines = f.readlines()
    struct_align = lines[1].rstrip()
    lib_id = lines[2].rstrip().split('>')[-1].split('.')[0]
    lib_align = lines[3].rstrip()
    return struct_align, lib_id, lib_align


def main(fasta_alignment_pair, pdb_file, dna_start_position, dna_chain, f5, f3):
    structure_align, lib_id, lib_align = parse_paired_alignment(fasta_alignment_pair)

    if '-' in structure_align:
        raise Exception('there is a gap in the structure sequence, which is not allowed. library sequence: '
                        '{0} was skipped\n'.format(lib_id))
    leading_gap_len = len(lib_align) - len(lib_align.lstrip('-'))
    trailing_gap_len = len(lib_align) - len(lib_align.rstrip('-'))
    if leading_gap_len > len(f5) or trailing_gap_len > len(f3):
        raise Exception('the leading or trailing gap in {0} exceeds the '
                        'length of the given flanking dna sequences\n'.format(lib_id))
    if '-' in lib_align.lstrip('-').rstrip('-'):
        raise Exception('there is an internal gap in {0}, which is not allowed\n'.format(lib_id))

    if leading_gap_len > 0:
        lib_seq = f5[-leading_gap_len:] + lib_align[leading_gap_len:]
    if trailing_gap_len > 0:
        lib_seq = lib_seq[:-trailing_gap_len] + f3[:trailing_gap_len]

    assert len(lib_seq) == len(structure_align), 'adjusted library seq {0} is not the same length as ' \
                                                 'the structure seq\n'.format(lib_id)

    dna_abbreviation_dict = {'A': 'ADE', 'C': 'CYT', 'G': 'GUA', 'T': 'THY'}
    dna_def_list = []
    for index, base_tup in enumerate(zip(structure_align, lib_seq)):
        if base_tup[0] != base_tup[1]:
            position = dna_start_position + index
            dna_def_str = '{0}.{1}.{2}'.format(dna_chain, position, dna_abbreviation_dict[base_tup[1]])
            dna_def_list.append(dna_def_str)
    log_file = '{0}_{1}.log'.format(pdb_file.split('.')[0], lib_id)
    command = '/scratch/nwhoppe/rosetta_src_2015.38.58158_bundle/main/source/bin/rosettaDNA.linuxgccrelease ' \
              '-database /scratch/nwhoppe/rosetta_src_2015.38.58158_bundle/main/database/ ' \
              '@resfile.flags -in::file::s {0} -out:suffix _{1} ' \
              '-parser::protocol pr.xml -parser::script_vars ' \
              'dna_seq={2} > {3}\n'.format(pdb_file, lib_id, ','.join(dna_def_list), log_file)
    sys.stdout.write(command)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script which writes a shell script for running rosetta dna
    binding energy calculations. input includes: pdb structure with protein and dna, fasta file with sequence of dna in
    pdb and the flanking dna sequence flanking the sequence found in the pdb,
    fasta file with sequences for which binding energy calculation is desired. Flanking dna sequence is used to fill
    in terminal gaps in the alignment between a query dna sequence and the bound dna sequence. no internal gaps are
    allowed. for now, both an aligned replacement and an unaligned replacement will be made""")

    parser.add_argument('-p', '--pdb', default='3lsr_native.pdb', help='name of pdb structure used as input to rosetta')
    parser.add_argument('-s', '--start_pos', type=int, default=5, help='pdb number of first dna base')
    parser.add_argument('-c', '--dna_chain', default='C', help='chain of top strand of dna in pdb')
    parser.add_argument('-f5', default='TTGACA', help='bases 5 prime of library sequences')
    parser.add_argument('-f3', default='TATAAT', help='bases 3 prime of library sequences')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-a', '--alignment', required=True,
                          help='two aligned sequences in fasta format. file should be four lines. first two lines are'
                               ' the dna sequence found in the pdb. second two lines are library sequence that is to '
                               'be modeled into the pdb with Rosetta')

    args = parser.parse_args()
    main(args.alignment, args.pdb, args.start_pos, args.dna_chain, args.f5, args.f3)
