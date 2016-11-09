# python2 script for converting protein to DNA sequence and checking for restriction enzyme digest sites

from Bio import SeqIO
from Bio import Restriction
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from CodonUsage import sorted_codon_table
import re
import sys

base_file = sys.argv[2]
new_file = base_file[:-6] + "_fixedseqs.fasta"
print new_file
oligo_file = base_file[:-6] + "_oligos.fasta"


# function for reading out the fasta sequences from the fasta file.  It converts the file to a list
def open_fastas(fasta_file, sequence_type):
    with open(fasta_file, 'rU') as the_file:
        sequences = SeqIO.parse(the_file, 'fasta', sequence_type)
        final_seqs = []
        for sequence in sequences:
            final_seqs.append(sequence)
    return final_seqs


# -----------------------CHECKS---------------------------------------
# function to make sure the wt dna sequence encodes the right protein and has a stop codon
def check_native_dna(dna_seq):
    if dna_seq.seq.translate().find('*') == -1:  # the .find() function returns -1 if the search is not found
        print "Stop codon not found..."
        return False
    elif (len(dna_seq.seq)) / 3 != (dna_seq.seq.translate().find('*') + 1):
        # accounts for python counting for indices but not for lengths
        # length of wt dna sequence divided by 3 = number of amino acids in full length protein; the find function returns lowest index containing the search function
        print "Premature stop codon"
        return False
    else:
        return True


# returns True/False depending on if a protein sequence has a stop codon
def stop_codon_check(protein_seq):
    if not protein_seq.seq.endswith('*'):
        #		print protein_seq
        #		print "No stop codon"
        return False
    else:
        return True


# for adding stop codon by checking through a protein sequence
def add_stop(protein_seq):
    #	print "add_stop check"
    #	print protein_seqs
    if stop_codon_check(protein_seq) == False:
        print "no stop codon"
        protein_seq = protein_seq + Seq('*')
    else:
        pass
    return protein_seq


def length_check(wt_seq, protein_seq):
    if len(wt_seq.seq.translate()) != len(seq.seq):
        print "Lengths not the same"
        return False
    else:
        return True


# input is a single sequence, checks to make sure only natural amino acids foud
def aa_check(protein_seq):
    aa_list = 'ACDEFGHIKLMNPQRSTVWY*'
    for aa in protein_seq.seq:
        if aa not in aa_list:
            print aa
            #			print protein_seq.seq[aa]
            return False
        else:
            return True


def dna_to_protein_check(protein_dna, protein_aas):
    translation = protein_dna.seq.translate()
    if str(translation) != str(protein_aas.seq):
        return False
        print "Protein DNA translation and amino acid sequence don't match"
    else:
        return True


# ------------------------------CHECKS END------------------------------

# gives as a dictionary of enzyme: cut position starts given a DNA seqobj
def check_res_site(DNA_seq):
    rb = RestrictionBatch([BsaI, EcoRI, NheI])
    print "Restriction site searched: ", rb
    cut_locations = rb.search(DNA_seq.seq)
    print "cut_locations: ", cut_locations
    return cut_locations


# gives a list of tuples of the (start, end) locations of the cut sites given a dna seqobj and the dictionary of enzyme: cut locations (from check_res_site())
def find_cut_pos(mutated_dna_seqobj, cut_loc):
    cut_startend = []
    for enzyme in cut_loc:
        if cut_loc[enzyme]:  # exists...
            if enzyme == BsaI:
                restriction_enz = enzyme
                fwd_seq = str(mutated_dna_seqobj.seq.upper())
                rev_seq = str(mutated_dna_seqobj.seq.reverse_complement().upper())
                #				print "5'-", fwd_seq, "-3'", "--fwd strand"
                #				print "5'-", rev_seq, "-3'", "--rev strand"
                #				print restriction_enz.site
                for position in re.finditer(str(restriction_enz.site), fwd_seq):
                    start_nt = position.start()
                    end_nt = position.end()
                    print "Cut site found on the forward strand between positions ", start_nt, " and ", end_nt
                    cut_startend.append((start_nt, end_nt))
                # appends the nt start # and end # as a tuple to the cut_startend list for the forward sequence
                for position in re.finditer(str(restriction_enz.site), rev_seq):
                    start_nt = len(
                        mutated_dna_seqobj.seq) - position.start()  # counts from same end as + strand, but the start nt would be > end nt
                    end_nt = len(mutated_dna_seqobj.seq) - position.end()  # see above
                    print "Cut site found on reverse strand between positions ", start_nt, " and ", end_nt
                    cut_startend.append((start_nt, end_nt))
            else:
                restriction_enz = enzyme
                fwd_seq = str(mutated_dna_seqobj.seq.upper())
                #				print "5'-", fwd_seq, "-3'", "--fwd strand"
                #				print restriction_enz.site
                for position in re.finditer(str(restriction_enz.site), fwd_seq):
                    start_nt = position.start()
                    end_nt = position.end()
                    print "Cut site found on the forward strand between positions ", start_nt, " and ", end_nt
                    cut_startend.append((start_nt, end_nt))
                #			print cut_startend

    return cut_startend


# generates a list of synonymous codons given the input of a codon
def alternative_codons(old_codon):
    other_codons = []
    for AA, codons in sorted_codon_table.items():
        if str(old_codon) in codons:
            for item in sorted_codon_table[AA]:
                other_codons.append(item)
            #			print "#######", sorted_codon_table[AA]
            other_codons.remove(old_codon)  # removes old codon from the list
    return other_codons


# function to fix restriction sites.  the inputs are the native_protein dna sequence, the mutated dna sequence obj, and the list of cut positions given by the find_cut_pos function
# Note: it looks like native_prot not required
def fix_restriction(native_prot, mutated_seq, cut_pos):
    for items in cut_pos:
        test = mutated_seq.seq.upper()
        other_codons = []
        #		print items
        start_nt = items[0]  # first item in tuple
        end_nt = items[1]  # second item in tuple
        quotient, remainder = divmod(start_nt, 3)[0], divmod(end_nt, 3)[
            1]  # divmod gives the quotient remainder as 2 parts of a list, assigns each to the proper variable
        codon_start = quotient * 3  # gives earliest codon containing the cut site
        if start_nt < end_nt:
            # the cut site is on the + strand
            print "Cut site on + strand"
            if codon_start <= start_nt:
                # The start of the restriction site is past the first codon containing it
                codon_start += 3
                # adjust so the codon changed is in the middle of the restriction site
            codon_end = codon_start + 3
            # THIS PART NECESSARY? most cut sites are palindromic--I wonder if this screwed up the old script trying to fix the same site on the + and the - strand
            # ONLY necessary for the BsaI (GGTCTC)
        elif start_nt > end_nt:
            # cut site is on - strand
            print "Cut site on - strand"
            if codon_start >= start_nt:
                codon_start = codon_start - 3
            codon_end = codon_start + 3
        else:
            print "start = end, exiting"
            exit()
            # adjusting if codon is the last one in the sequence (which it shouldn't be because that would be the stop codon...
        if len(test) == codon_end:
            codon_start = codon_start - 3
            codon_end = codon_end - 3

        codon_to_replace = test[codon_start:codon_end]
        #		print codon_to_replace
        new_codon = "TAA"  # stop codon as default replacement if no other is found
        other_codons = alternative_codons(codon_to_replace)
        #		print other_codons
        while True:
            if len(other_codons) == 0:
                print "Error: No alternative codons found, changing codon to replace"
                if len(test) - 3 == codon_end:
                    codon_start = codon_start - 3
                    codon_end = codon_end - 3
                    codon_to_replace = test[codon_start:codon_end]
                    other_codons = alternative_codons(codon_to_replace)
                else:
                    codon_start = codon_start + 3   ### this could give a codon past the restriction site if codon_start = start_nt
                    # and the restriction site is not longer than 6 nucleotides
                    codon_end = codon_end + 3
                    codon_to_replace = test[codon_start:codon_end]
                    other_codons = alternative_codons(codon_to_replace)

            else:
                for codon in other_codons:
                    new_codon = str(codon)
                    if str(codon_to_replace) != new_codon:
                        break
                break
        print "REMARK: replacing old codon ", codon_to_replace, "with synonymous new codon ", new_codon
        upstream_seq = test[:codon_start]
        #		print "Upstream_seq: ", upstream_seq
        downstream_seq = test[codon_end:]
        #		print "Dnstream_seq: ", downstream_seq
        temp = upstream_seq + new_codon + downstream_seq
        #		print "Temp: ", temp
        mutated_seq = SeqRecord(Seq("".join(temp), IUPAC.unambiguous_dna), id=mutated_seq.id, name=mutated_seq.name,
                                description="codon optimized for mutations, restriction sites checked")
    return mutated_seq


def convert_to_DNA(native_prot_dict, mutated_seq_obj):
    variant_dna_codons = []
    for i in range(len(mutated_seq_obj.seq)):
        codon = str(native_prot_dict[i + 1][
                        1])  # dictionary starts at 1, not 0 (i+1), and the second [1] indicates the dna seq, not amino acid letter
        #		print codon
        nat_aa = str(native_prot_dict[i + 1][0])
        if mutated_seq_obj.seq[i] != nat_aa:
            if mutated_seq_obj[i] == '-':
                variant_dna_codons.append(codon)
            else:
                pass
            codon = sorted_codon_table[str(mutated_seq_obj.seq[i])][0]
        #			print codon
        #			print sorted_codon_table[str(mutated_seq_obj.seq[i])]
        variant_dna_codons.append(str(codon))
    variant_dna = "".join(variant_dna_codons)
    #	print variant_dna
    new_mut_seq = Seq(variant_dna, IUPAC.unambiguous_dna)
    new_mut_dna_seq_obj = SeqRecord(new_mut_seq, id=mutated_seq_obj.id, name=mutated_seq_obj.name,
                                    description=mutated_seq_obj.description)
    #	print new_mut_seq_obj
    #	cut_locations = {}
    cut_locations = check_res_site(new_mut_dna_seq_obj)
    #	print cut_locations
    FoundResSite = False
    for item in cut_locations:
        if cut_locations[item]:  # exists...
            print cut_locations[item]
            FoundResSite = True
            print "WARNING:", item, "cut site found in sequence ", new_mut_dna_seq_obj.name
    else:
        pass
    if FoundResSite == True:  # is true...
        cut_pos = find_cut_pos(new_mut_dna_seq_obj, cut_locations)
        fixed_mut_dnaseqobj = fix_restriction(native_dna_seq.seq, new_mut_dna_seq_obj, cut_pos)
    #		print fixed_mut_dnaseqobj
    #	print "aaaaaaa", cut_locations.values()
    #	print cut_locations
    # loop to account for changes introducing more cut sites
    while cut_locations.values() != [] and FoundResSite == True:
        cut_locations = check_res_site(fixed_mut_dnaseqobj)
        FoundResSite = False
        for item in cut_locations:
            if cut_locations[item]:  # exists...
                print cut_locations[item]
                FoundResSite = True
                print "WARNING:", item, "cut site found in sequence ", new_mut_dna_seq_obj.name
        else:
            #        		print "No restriction site found for", item
            pass
        if FoundResSite == True:  # is true...
            cut_pos = find_cut_pos(fixed_mut_dnaseqobj, cut_locations)
            fixed_mut_dnaseqobj = fix_restriction(native_dna_seq.seq, fixed_mut_dnaseqobj, cut_pos)
            print cut_locations
        else:
            print "No restriction sites found"
            break
    if FoundResSite == True and not dna_to_protein_check(fixed_mut_dnaseqobj, mutated_seq_obj):
        print "Error: translated fixed sequence does not match original protein"
        print exit()
    else:
        try:
            dna_to_protein_check(fixed_mut_dnaseqobj, mutated_seq_obj)
            # will give unbound local error if fixed_mut_dnaseqobj is not defined (which it isn't if there are no changes to the sequence via the fixrestrictionsite function)
        except UnboundLocalError:
            print "Returning Old Sequence"
            return new_mut_dna_seq_obj

        else:
            print "Returning Fixed Sequence"
            return fixed_mut_dnaseqobj


# added to account for missing M, input is the list of seqobj and the list containing the wt protein seqobj
def seq_addition_tetr(record, wt_prot_list):
    wt_prot = []
    for prot in wt_prot_list:
        for aa in prot.seq.translate():
            wt_prot.append(aa)
    changed_record = []
    for sequence in record:

        #		print(len(sequence.seq))
        new_sequence = ''
        for letter in sequence.seq:
            if letter == "-":
                the_index = len(new_sequence)
                new_sequence = new_sequence + wt_prot[the_index]
            else:
                new_sequence = new_sequence + letter
        if '4AC0' in sequence.name:
            new_sequence = new_sequence + 'ESGS*'
        else:
            pass
        #		print new_sequence
        #		print len(sequence.seq)
        testing = Seq(new_sequence, IUPAC.protein)
        test_record = SeqRecord(testing, id=sequence.id, name=sequence.name, description=sequence.description)
        changed_record.append(test_record)
    return changed_record


# function for getting the sequence fragment of what we're looking for (B0 or B1)
# just one sequence
# primers taken from megan's text document "flankingseqs_MphR.txt"
def get_fragments(sequence):
    BsaI_5 = 'GGTCTC'
    BsaI_3 = 'GAGACC'

    if '4AC0' in sequence.name:
        if 'B0' in sequence.name:
            print "4ac0_b0"
            #			F_prim_add = 'CCCGTCCCTG'
            F_prim_add = 'gtgacccgtccctg'  # removed 5'nt to make total length 170
            #			R_prim_add = 'gggcagaggt'
            R_prim_add = 'gggcagaggtcgac'  # removed 3'nt to make total length 170
            five_Flank = 'AAGAT'
            three_Flank = 'gcctt'
            start_nt = 77 * 3  # starts at F78, but starts count from 0
            # note, need to make sure these are the right sequences....
            end_nt = 117 * 3  # ends at 117 (not inclusive)
            #			print start_nt
            #			print end_nt
            fragment = sequence.seq[start_nt:end_nt]
            #			print fragment
            #			print fragment.translate()
            #			print oligo


            """
            F_primer_seq = Seq(F_prim_add+BsaI_5, IUPAC.unambiguous_dna)
            R_primer = Seq(R_prim_add+BsaI_3, IUPAC.unambiguous_dna)
            R_primer_seq = R_primer.reverse_complement()
            print "forward primer: ", F_primer_seq
            print R_primer_seq
            """

        elif 'B1' in sequence.name:
            print "4ac0_b1"
            #			F_prim_add = 'tgcccgctgt'
            F_prim_add = 'tgcccgctgtcttca'
            #			R_prim_add = 'cccggcagtg'
            R_prim_add = 'tgtagcccggcagtg'
            five_Flank = 'AAGTA'
            three_Flank = 'CATTT'
            start_nt = 99 * 3
            end_nt = 138 * 3
            fragment = sequence.seq[start_nt:end_nt]
        #			print fragment
        #			print fragment.translate()
        #			print oligo

        else:
            print sequence.name
            print "Sequence not recognized"
    elif "2uxo" in sequence.name:
        if 'B0' in sequence.name:
            print "2uxo_B0"
            F_prim_add = 'CGATCGTGCCCACCT'
            R_prim_add = 'AGTTGGAGCCCGCAC'
            five_Flank = 'cactg'
            three_Flank = 'gttct'
            start_nt = 62 * 3  # starts at 63, but starts count from 0
            # note, need to make sure these are the right sequences....
            end_nt = 100 * 3  # ends at 100 (not inclusive)
            #			print start_nt
            #			print end_nt
            fragment = sequence.seq[start_nt:end_nt]
        #			print fragment
        #			print fragment.translate()
        #			print oligo
        elif 'B1' in sequence.name:
            print "2uxo_B1"
            F_prim_add = 'CTGGTGCGTCGTCT'  # remove 1nt from 5' end
            R_prim_add = 'GGCGAACACTTCCC'  # remove 1nt from 3' end
            five_Flank = 'tggat'
            three_Flank = 'cgttg'
            start_nt = 136 * 3  # starts at F78, but starts count from 0
            # note, need to make sure these are the right sequences....
            end_nt = 176 * 3  # ends at 117 (not inclusive)
            #			print start_nt
            #			print end_nt
            fragment = sequence.seq[start_nt:end_nt]
        #			print fragment
        #			print fragment.translate()
        #			print oligo
        else:
            print sequence.name
            print "Sequence not recognized"
    oligo = F_prim_add + BsaI_5 + five_Flank + fragment + three_Flank + BsaI_3 + R_prim_add
    return oligo


"""
test = open_fastas(sys.argv[1], IUPAC.unambiguous_dna)
print test
for seq in test:
	lel = check_res_site(seq)
	the_list = find_cut_pos(seq, lel)
	
	fix_restriction("what", seq, the_list)
"""

wt = open_fastas(sys.argv[1], IUPAC.unambiguous_dna)
assert len(wt) == 1, 'Only one wildtype fasta can be specified and it must be a DNA sequence'

for item in wt:
    native_dna_seq = item

if not check_native_dna(native_dna_seq):
    print "Make sure the wt seq has a stop codon"
    exit()

native_prot = {}
for i in range(3, len(native_dna_seq) + 1, 3):
    if not (i / 3) in native_prot:
        native_prot[i / 3] = ()
    native_prot[i / 3] = (native_dna_seq.seq[i - 3:i].translate(), native_dna_seq.seq[i - 3:i])
# print native_prot
# native_prot is a dictionary of {#: protein aa, codon} as seqobj

mut = open_fastas(sys.argv[2], IUPAC.protein)
mut = seq_addition_tetr(mut, wt)

mutated_prot_seqs = []

# ---------------CHECK IMPLEMENTATION----------------



for seq in mut:
    #	print seq
    seq = add_stop(seq)
    if not length_check(native_dna_seq, seq):
        print "Protein lengths do not match for sequence", seq.id
        exit()
    #	print seq.seq
    if not aa_check(seq):
        print "Unnatural amino acids found for sequence", seq.id
        exit()
    mutated_prot_seqs.append(seq)

# -----------------OBTAINING FULL SEQUENCES--------------------------

oligo_seqs = []
for seq in mutated_prot_seqs:
    dna_output = convert_to_DNA(native_prot, seq)
    print dna_output

    with open(new_file, 'a') as output:
        output.write(">%s %s\n" % (dna_output.name, dna_output.description))
        output.write("%s\n" % dna_output.seq)
    oligos = get_fragments(dna_output)
    oligo_seqs.append(oligos)
    with open(oligo_file, 'a') as output:
        output.write(">%s\n" % (dna_output.name))
        output.write("%s\n" % oligos)
    print "----------"
print "-----------"
for item in oligo_seqs:
    print item
    print len(item)
# NOTE: NEED TO FIX THE APPEND SYSTEM TO NOT REWRITE THE SCRIPT EVERY TIME.  its because i modified the script to run on just a seq, need a global list in the for string.
