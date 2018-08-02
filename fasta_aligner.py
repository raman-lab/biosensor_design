# script using biopython to pull out fastas, align them in a pssm, and pull out the residues that were changed compared to the wild type
#knishikawa edit 20160606 for addition of mutation_counter function
from sys import argv
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import os
import re
import matplotlib.pyplot as plt
import numpy as np
script, wild_type_fasta = argv
#*_modified.fasta

#print aligned_sequences
#print aligned_sequences[1] #pulls information on a particular sequence
#print aligned_sequences[1].seq #can get just the sequence
#print aligned_sequences[1].name #pulls sequence name
#print aligned_sequences[1][0] #just residue
#initial_align = AlignInfo.SummaryInfo(aligned_sequences)
#print initial_align

# converts the fasta sequences to dashes when the sequences match that of the wild-type sequence
def dashed_alignment(multiple_seq_alignment):
	wild_file = open(wild_type_fasta, 'rU') #opens the wild type fasta sequence given as wild_type_fasta in argv
	wt_data = SeqIO.read(wild_file, 'fasta') #throws the sequence into a SeqRecord object
	#wt_seq = wt_data.seq #may not be necessary
	new_fasta = open(dashed_fasta, 'w') #opens the third argv file to place the dashed sequences
	
	for sequence in multiple_seq_alignment:
		new_sequence = "" #creates a new list to hold the sequence
		sequence_name = sequence.name #gives the name of the sequence in the fasta file a variable
		#stripped_seq_name = strip_punctuation(sequence_name)
		mod_sequence_name ='>'+sequence_name+'\n' #changes the sequence name into fasta format
		new_fasta.write(mod_sequence_name) #writes the new fasta-compatible sequence name to a file without quotations
		for res_num in range(0, len(sequence)): #passes through each residue in the sequence by number
			residue = sequence[res_num] #sets the amino acid at each position to a variable
			mod_residue = residue.strip('"\'') #strips quotations from the amino acid (normally found as 'A'), probably not necessary given the write function later
			if wt_data[res_num] == residue: #matches the wt sequence at position [res_num] to the residue currently being examined
				new_sequence = new_sequence+'-' #places a dash in the new empty sequence if they match
			else:
				new_sequence = new_sequence+mod_residue #if the amino acid at this position does not match the wt sequence, the residue is placed in there instead
		mod_sequence = new_sequence+'\n' #places a newline character after the completed sequence is created
			
		new_fasta.write(mod_sequence) #writes the mod sequence such that the sequence is inserted without quotations that mess with fasta files
	wild_file.close()
	new_fasta.close()
	#return new_fasta




# this function generates a pssm for the dashed sequences completed in an earlier function
def pssm_generator():
	wild_file = open(wild_type_fasta, 'rU') #opens the wild type fasta sequence given as wild_type_fasta in argv
	wt_data = SeqIO.read(wild_file, 'fasta') #throws the sequence into a SeqRecord object
	dashed_file = open(dashed_fasta, 'rU') #opens the third argv variable that was written in the dashed_alignment function
	dashed_aligned_sequences = AlignIO.read(dashed_file, 'fasta') #throws the dashed alignment into Biopython
	summary_align = AlignInfo.SummaryInfo(dashed_aligned_sequences) #required for pssm generation
	#print summary_align
	consensus_sequence =  summary_align.dumb_consensus() #gets consensus sequence from the summary align information
	consensus_sequence = str(consensus_sequence)
	print consensus_sequence
	consensus_file = open(new_consensus, 'w')
	consensus_name = '>consensus_' + wt_data.name + '\n'
	consensus_file.write(consensus_name)
	consensus_file.write(consensus_sequence)
	#print consensus_sequence
	my_pssm = summary_align.pos_specific_score_matrix() #creates pssm from summary align information
	print my_pssm
	dashed_file.close()
	consensus_file.close()
	return my_pssm




#this function creates the graphs for amino acid abundance at each position that differ from teh wildtype in any sequence in the original list
def generate_graphs(my_pssm):
	dashed_file = open(dashed_fasta, 'rU')
	dashed_aligned_sequences = AlignIO.read(dashed_file, 'fasta') #this line and the previous line are required for getting the total number of sequences for later
	residue_number=0 #so we can iterate through positions
	for position in my_pssm:
		residue_number = residue_number+2 #converts 'programming counting' to the residue number (instead of starting at position 0 by count, will start at position 1)
		if position['-'] != len(dashed_aligned_sequences): #this is so that we don't look at positions were all the sequences were a dash (same as the wild-type)
		
			print residue_number
			print position.keys() #prints amino acids
			print position.values() #prints abundance of each amino acids
			num_amino_acids = len(position.values()) #prints the number of amino acids
		
			fig, ax = plt.subplots()
		#graphing stuff
			graph = ax.bar(range(len(position.values())), position.values(), width=1, color='b') #bar(xvalues, yvalues, width, color), x values must be a number
			plt.xticks(range(len(position.values())), position.keys()) #converts the number on the xaxis to a letter
			ax.set_xlabel('Residue', fontsize = 20)
			ax.set_ylabel('Abundance', fontsize = 20)
			ax.set_title(residue_number, fontsize = 30)
			plt.savefig('TetR_Design_Position_%r.png' %residue_number)
		else:
			print residue_number
	dashed_file.close()
		

#generate_graphs(my_pssm)
def mutation_counter(multiple_seq_alignment):
	wild_file = open(wild_type_fasta, 'rU') #opens the wild type fasta sequence given as wild_type_fasta in argv
	wt_data = SeqIO.read(wild_file, 'fasta') #throws the sequence into a SeqRecord object
	#wt_seq = wt_data.seq #may not be necessary
	mutations_counter = []
	seq_counter = 0
	for sequence in multiple_seq_alignment:
		seq_counter += 1
		num_mut = 0
		for res_num in range(0, len(sequence)): #passes through each residue in the sequence by number
			residue = sequence[res_num] #sets the amino acid at each position to a variable
			mod_residue = residue.strip('"\'') #strips quotations from the amino acid (normally found as 'A'), probably not necessary given the write function later
			if wt_data[res_num] == residue: #matches the wt sequence at position [res_num] to the residue currently being examined
				pass
			else:
				num_mut += 1
				
		mutations_counter.append(num_mut)
		
	xaxisbin = np.arange(0, 12, 1)
	gtitle = "Mutation_count_" + base_file
	plt.ylim([0.0, seq_counter]) #changed to plt.ylim(lower,upper) for matplotlib 0.99
	plt.hist(mutations_counter, bins=xaxisbin)
	plt.xlabel("Num. Mutations")
	plt.ylabel("Abundance")
	plt.title(gtitle)
	#plt.axis([-300, -240])
	plt.savefig(gtitle + "_" + ".png")
	plt.clf()
		
		



currentpath = os.getcwd()
for file in os.listdir(currentpath):
	if file.endswith('_modified.fasta'):
		base_file = file[:-15]
		dashed_fasta=base_file+'_dashed.fasta'
		new_consensus=base_file+'_consensus.fasta'

		the_file = open(file, 'rU') #opens first file listed in argv
		fasta_sequences = SeqIO.parse(the_file, 'fasta') #pulls out sequences for use
		aligned_sequences = AlignIO.read(the_file, 'fasta') #aligns sequences from the original fasta file
		print "Adjusting fasta alignment to dashes..."
		dashed_alignment(aligned_sequences)
		print "done"
		print "Aligning dashed sequences..."
		my_pssm = pssm_generator() #sets the pssm generated by this function to the variable my_pssm
		mutation_counter(aligned_sequences)
		the_file.close()

	else:
		pass


