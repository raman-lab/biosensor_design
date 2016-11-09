#!/usr/bin/env python

import argparse
import copy
import re
import sys
import collections
from Bio import SeqIO
from Bio import Restriction
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from CodonUsage import sorted_codon_table


def parse_wt_sequences(wt_dna_fasta):
    """load in wt sequence and return wt dna sequence and dictionary of the form
    {amino acid position: amino acid, codon}
    checks for: existence and position of stop codon"""
    with open(wt_dna_fasta, 'r') as f:
        sequence_generator_obj = SeqIO.parse(f, 'fasta', IUPAC.unambiguous_dna)
        sequences = []
        for sequence_obj in sequence_generator_obj:
            sequences.append(sequence_obj.seq.upper())
    assert len(sequences) == 1, 'Only one wild type dna sequence can be present in the wild type dna fasta file'
    wt_dna_seq = sequences[0]
    if wt_dna_seq.translate().find('*') == -1:
        raise Exception("Stop codon not found in wild type sequence")
    elif (len(wt_dna_seq)) / 3 != (wt_dna_seq.translate().find('*') + 1):
        raise Exception("Premature stop codon found in wild type sequence")

    native_protein_dict = {}
    for i in range(3, len(wt_dna_seq) + 1, 3):
        if not (i / 3) in native_protein_dict:
            native_protein_dict[i / 3] = ()
        native_protein_dict[i / 3] = (wt_dna_seq[i - 3:i].translate(), wt_dna_seq[i - 3:i])
    return wt_dna_seq, native_protein_dict


def parse_protein_sequences(protein_seq_files):
    """load in protein amino acid sequences from fasta files"""
    all_sequences_objs = []
    for fasta in protein_seq_files:
        with open(fasta, 'r') as f:
            sequences_in_fasta = SeqIO.parse(f, 'fasta', IUPAC.protein)
            for sequence_obj in sequences_in_fasta:
                all_sequences_objs.append(sequence_obj)
    return all_sequences_objs


def convert_to_dna(protein_sequence, wt_protein_dict):
    """converts protein variant sequence to dna
    if a dash is present in the variant sequence, it is replaced by the wt amino acid"""
    variant_dna_codons = []
    for index in range(0, len(protein_sequence.seq)):
        wt_aa = str(wt_protein_dict[index + 1][0])
        codon = str(wt_protein_dict[index + 1][1])
        variant_aa = protein_sequence.seq[index]
        if variant_aa != wt_aa:
            if variant_aa is not '-':
                codon = sorted_codon_table[str(variant_aa)][0]
        variant_dna_codons.append(codon)
    variant_dna_str = "".join(variant_dna_codons)
    variant_dna_seq = Seq(variant_dna_str, IUPAC.unambiguous_dna)
    variant_dna_seq_obj = SeqRecord(variant_dna_seq, id=protein_sequence.id, name=protein_sequence.name,
                                    description=protein_sequence.description)
    return variant_dna_seq_obj


def add_flanking_nucleotides(dna_variant, protein_variant, wt_dna, primer_file, amino_acid_range, restriction_enzyme):
    """from a variant dna sequence, the following operations are performed:
        desired range is extracted, 5 wild type nucleotides are added on both ends,
        restriction enzyme site is added, specified primers are added
    depending on the length of the resulting oligo, nucleotides will be added or will be removed
    from the 5 prime end of the forward primer to make len(sequence) % 3 == 2 true"""

    # get primers to be added later. primers in the else statement are known to work tetR and ttgR oligos and will
    # be used in order of tetR B0, tetR B1, ttgR B0, ttgR B1
    if primer_file:
        five_primers = []
        three_primers = []
        with open(primer_file, 'r') as f:
            for line in f:
                split_line = line.split()
                column1_primer = Seq(split_line[0], IUPAC.unambiguous_dna)
                column2_primer = Seq(split_line[1], IUPAC.unambiguous_dna)
                reverse_complement_reverse_primer = column2_primer.reverse_complement()
                five_primers.append(column1_primer)
                three_primers.append(reverse_complement_reverse_primer)
    else:
        five_primers = [
            Seq("GGGTCACGCGTAGGA", IUPAC.unambiguous_dna),
            Seq("TGCCCGCTGTCTTCA", IUPAC.unambiguous_dna),
            Seq("CGATCGTGCCCACCT", IUPAC.unambiguous_dna),
            Seq("ACTGGTGCGTCGTCT", IUPAC.unambiguous_dna),
        ]
        three_primers = [
            Seq("GTGTGGCTGCGGAAC", IUPAC.unambiguous_dna),
            Seq("TGTAGCCCGGCAGTG", IUPAC.unambiguous_dna),
            Seq("AGTTGGAGCCCGCAC", IUPAC.unambiguous_dna),
            Seq("GGCGAACACTTCCCG", IUPAC.unambiguous_dna),
        ]

    # get dna_fragment from desired range
    if amino_acid_range:
        amino_acids = amino_acid_range.split('-')
        first_codon = int((amino_acids[0] - 1) * 3)
        last_codon = int(amino_acids[1] * 3)
    else:
        name = dna_variant.name.upper()
        if '4AC0' in name and 'B0' in name:
            first_codon = 77 * 3
            last_codon = 117 * 3
            five_primer = five_primers[0]
            three_primer = three_primers[0]
        elif '4AC0' in name and 'B1' in name:
            first_codon = 99 * 3
            last_codon = 138 * 3
            five_primer = five_primers[1]
            three_primer = three_primers[1]
        elif '2UXO' in name and 'B0' in name:
            first_codon = 62 * 3
            last_codon = 100 * 3
            five_primer = five_primers[2]
            three_primer = three_primers[2]
        elif '2UXO' in name and 'B1' in name:
            first_codon = 136 * 3
            last_codon = 176 * 3
            five_primer = five_primers[3]
            three_primer = three_primers[3]
        else:
            raise Exception('if amino acid range is not specified, 4AC0 or 2uxo and B0 or B1 must be in '
                            'each sequence name')
    protein_fragment = protein_variant.seq[first_codon / 3: last_codon / 3]
    dna_fragment = dna_variant.seq[first_codon:last_codon]

    # add 5 wild type nucleotides
    dna_fragment_w_buffer = Seq('{0}{1}{2}'.format(wt_dna[first_codon - 5:first_codon], dna_fragment,
                                                   wt_dna[last_codon:last_codon + 5]), IUPAC.unambiguous_dna)

    # add restriction sites
    restriction_batch = Restriction.RestrictionBatch([restriction_enzyme])
    five_restriction_site = Seq(restriction_batch.get(restriction_enzyme).site, IUPAC.unambiguous_dna)
    three_restiction_site = five_restriction_site.reverse_complement()
    dna_fragment_wo_primer = Seq('{0}{1}{2}'.format(five_restriction_site, dna_fragment_w_buffer,
                                                    three_restiction_site), IUPAC.unambiguous_dna)

    # add primers
    oligo = Seq('{0}{1}{2}'.format(five_primer, dna_fragment_wo_primer, three_primer), IUPAC.unambiguous_dna)

    # manipulate oligo to put it in frame with Kan in construct and make it compatible with 170 nucleotide limit
    assert len(oligo) < 173, 'Error: oligo must be 172 nucleotides or less. currently the length is {0}. this ' \
                             'is not possible if primers are 15 nucleotides, restriction sites are 6 nucleotides, ' \
                             'the buffer is 5 nucleotides, and the variable dna fragment is 120 nucletides. please ' \
                             'ensure these limits are adhered to'.format(len(oligo))
    mutable_oligo = oligo.tomutable()
    if len(oligo) == 172:
        mutable_oligo = mutable_oligo[1:-1]
    oligo = mutable_oligo.toseq()

    dna_fragment_seq_obj = SeqRecord(dna_fragment, id=dna_variant.id, name=dna_variant.name,
                                     description=dna_variant.description)
    oligo_seq_obj = SeqRecord(oligo, id=dna_variant.id, name=dna_variant.name, description=dna_variant.description)
    return oligo_seq_obj, dna_fragment_seq_obj, protein_fragment


def find_replace_cut_sites(restriction_sites, input_oligo, dna_fragment, restriction_enzyme, oligo_name):
    locations_to_fix = []
    for enzyme in restriction_sites:
        for site_type in restriction_sites[enzyme]:
            site_to_check = restriction_sites[enzyme][site_type]
            if str(enzyme) == restriction_enzyme:
                for occurrence in re.finditer(str(site_to_check), str(input_oligo)):
                    if occurrence.start() not in [16, len(input_oligo) - 21]:
                        locations_to_fix.append((occurrence.start(), occurrence.end()))
                        sys.stderr.write('Warning: restriction {0} for {1} found in oligo {2} at position {3}: {4}\n'.
                                         format(site_type, enzyme, oligo_name, occurrence.start() + 1, input_oligo))
            else:
                for occurrence in re.finditer(str(site_to_check), str(input_oligo)):
                    locations_to_fix.append((occurrence.start(), occurrence.end()))
                    sys.stderr.write('Warning: restriction {0} for {1} found in oligo {2} at position {3}: {4}\n'.
                                     format(site_type, enzyme, oligo_name, occurrence.start() + 1, input_oligo))
    mutable_oligo = input_oligo.tomutable()
    for location_to_fix in locations_to_fix:
        first_nucleotide = location_to_fix[0]
        first_codon = (first_nucleotide // 3) * 3
        if first_codon == first_nucleotide:
            first_codon += 3
            codon_to_replace = input_oligo[first_codon:first_codon + 3]
            alternate_codons_to_replace = [input_oligo[first_codon - 3: first_codon]]
        else:
            first_codon += 3
            codon_to_replace = input_oligo[first_codon:first_codon + 3]
            alternate_codons_to_replace = [input_oligo[first_codon - 3: first_codon],
                                           input_oligo[first_codon + 3: first_codon + 6]]
        if codon_to_replace in ['ATG', 'TGG']:
            # codon_to_replace corresponds to methionine or tryptophan, which only have one codon each
            # instead use alternate codon(s) found in the restriction site
            fail_count = 0
            for a, alternate in enumerate(alternate_codons_to_replace):
                if alternate in ['ATG', 'TGG']:
                    fail_count += 1
                else:
                    codon_to_replace = alternate
                    if a == 0:
                        codon_start = first_codon - 3
                    else:
                        codon_start = first_codon + 3
                    break
            assert fail_count != len(alternate_codons_to_replace), 'Fatal Error: all codons in the restriction site ' \
                                                                   'correspond to methionine or tryptophan, which ' \
                                                                   'only have one codon each. the restriction site ' \
                                                                   'cannot be corrected. See oligo {} postion {} ' \
                                                                   'to ensure this is correct'.format(oligo_name,
                                                                                                      first_nucleotide)
        else:
            codon_start = first_codon

        codon_table = copy.deepcopy(sorted_codon_table)
        amino_acid_to_replace = str(Seq(str(codon_to_replace), IUPAC.unambiguous_dna).translate())
        for codon in codon_table[amino_acid_to_replace]:
            if codon != codon_to_replace:
                new_codon = codon
                codon_table[amino_acid_to_replace].remove(new_codon)
                break
            else:
                codon_table[amino_acid_to_replace].remove(codon_to_replace)
        assert mutable_oligo[codon_start:codon_start + 3] == codon_to_replace, 'Error: the codon being replaced {0} ' \
                                                                               'does does not match the offending ' \
                                                                               'restriction site codon {1}. please ' \
                                                                               'inspect the oligo {2} at position ' \
                                                                               '{3}'.format(
            mutable_oligo[codon_start:codon_start + 3], codon_to_replace, input_oligo, codon_start)
        try:
            mutable_oligo[codon_start:codon_start + 3] = new_codon
        except NameError:
            raise Exception('A new codon has not been assigned from the codon table. please inspect the oligo {0} '
                            'at position {1}. the codon {2} at this position corresponds to amino acid {3}, '
                            'but was unable to be changed to an alternate codon from {4}\n'.format(
                input_oligo, codon_start, codon_to_replace, amino_acid_to_replace, codon_table[amino_acid_to_replace]))
    corrected_oligo = mutable_oligo.toseq()
    return corrected_oligo


def run_checks(oligo_seq_obj, dna_fragment_seq_obj, protein_fragment, restriction_enzyme):
    """oligo output checker. checks include:
     reading frame, stop condons in Kyle and Megan's backbone construction with Kan resistance as of 2016.07.19,
     restriction sites in oligo, equivalence of translated oligo fragment and variant amino acid block, uniqueness"""
    oligo = oligo_seq_obj.seq
    dna_fragment = dna_fragment_seq_obj.seq

    # restriction site check
    restriction_sites = collections.defaultdict(dict)
    restriction_batch = Restriction.RestrictionBatch(['BsaI', 'EcoRI', 'NheI', '{0}'.format(restriction_enzyme)])
    for enzyme in restriction_batch:
        site = Seq(restriction_batch.get(enzyme).site, IUPAC.unambiguous_dna)
        reverse_complement_site = site.reverse_complement()
        if site == reverse_complement_site:
            restriction_sites[enzyme] = {'site': site}
        else:
            restriction_sites[enzyme] = {'site': site, 'rc_site': reverse_complement_site}

    if len(oligo) == 170:
        input_oligo = Seq('CC{0}C'.format(str(oligo)), IUPAC.unambiguous_dna)
    else:
        input_oligo = Seq('C{0}'.format(str(oligo)), IUPAC.unambiguous_dna)

    found_cut_sites = True
    while found_cut_sites:
        fixed_oligo = find_replace_cut_sites(
            restriction_sites, input_oligo, dna_fragment, restriction_enzyme, oligo_seq_obj.name)
        if fixed_oligo == input_oligo:
            found_cut_sites = False
        else:
            input_oligo = fixed_oligo

    # reading frame and stop codons
    ecori_site = Restriction.EcoRI.site
    nhei_site = Restriction.NheI.site
    kan_construct = Seq('{0}{1}{2}A'.format(ecori_site, fixed_oligo, nhei_site), IUPAC.unambiguous_dna)
    assert len(kan_construct) % 3 == 0, 'Error: oligo causes kanamycin resistance gene in the assumed construct to ' \
                                        'be out of frame. check oligo {0}: {1}'.format(oligo_seq_obj.name,
                                                                                       str(fixed_oligo))
    translated_kan = kan_construct.translate()
    assert protein_fragment in translated_kan, 'Error: variable amino acid fragment {0} not found in the ' \
                                               'assumed construct {1}'.format(protein_fragment, translated_kan)
    assert translated_kan.find('*') == -1, 'Error: stop codon found in the reading frame of the assumed construct. ' \
                                           'dna sequence: {0} translated seq {1}'.format(str(kan_construct),
                                                                                         str(translated_kan))
    if len(fixed_oligo) == 173:
        temp_oligo = fixed_oligo[2:]
        final_oligo = temp_oligo[:-1]
    else:
        final_oligo = fixed_oligo[1:]
    return final_oligo


def make_oligos(protein_seq_files, wt_dna_fasta, amino_acid_range, primer_file, restriction_enzyme):
    """main function of script. used to control workflow"""
    wt_sequence, wt_protein_dict = parse_wt_sequences(wt_dna_fasta)
    protein_variants_objs = parse_protein_sequences(protein_seq_files)

    for variant in protein_variants_objs:
        dna_variant = convert_to_dna(variant, wt_protein_dict)

        oligo_seq_obj, dna_fragment_seq_obj, protein_fragment = add_flanking_nucleotides(
            dna_variant, variant, wt_sequence, primer_file, amino_acid_range, restriction_enzyme
        )
        checked_oligo = run_checks(oligo_seq_obj, dna_fragment_seq_obj, protein_fragment, restriction_enzyme)
        sys.stdout.write(">%s\n" % dna_variant.name)
        sys.stdout.write("%s\n" % checked_oligo)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""script to generate oligonucleotides from protein sequences
        script is highly specific to the construct used by Kyle and Megan
        run script with -h for list of options and default settings
        the structure of the output oligos will be:
            forward primer + cut site + nucleotide buffer + reverse transcribed DNA from amino acid range + nucleotide
            buffer + cut site + reverse complement of reverse primer"""
    )
    parser.add_argument("-r", "--restriction_enzyme", default='BsaI',
                        help="restriction enzyme that will be used to digest the oligo. "
                             "must be present in the Bio.Restriction library with a 6 nulceotide restriction site")
    parser.add_argument("-p", "--primer_file", help="file with list of primers to use in oligo generation. "
                                                    "primers at top of file will have priority. "
                                                    "if a primer introduces a stop codon, the next pair will be tried"
                                                    "file should have two columns: 1. forward primer 2. reverse primer "
                                                    "all primers should be 5'-3'"
                                                    "primers cannot contain the cut site of the input restriction"
                                                    "enzyme")
    parser.add_argument("-a", "--amino_acids",
                        help="inclusive range of amino acids from protein fasta files to be used in oligo generation. "
                             "i.e. 70-118 "
                             "defaults are programmed for the 'blocks' used in Rosetta designs protocols for tetR "
                             "and ttgR")
    requiredO = parser.add_argument_group('required arguments')
    requiredO.add_argument("-w", "--wild_type", help="wild type protein dna sequence in fasta format.")
    requiredO.add_argument("-f", "--fasta", nargs='*', required=True,
                           help="one or more protein fasta files from which amino acid range will be pulled")

    args = parser.parse_args()
    make_oligos(args.fasta, args.wild_type, args.amino_acids, args.primer_file, args.restriction_enzyme)
