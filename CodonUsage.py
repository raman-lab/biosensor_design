#!/usr/bin/env python

HighFreqCodonUsage = {
              'F':'TTT',
              'L':'CTG',
              'I':'ATT',
              'M':'ATG',
              'V':'GTT',
              'S':'AGC',
              'P':'CCG',
              'T':'ACC',
              'A':'GCG',
              'Y':'TAT',
              'H':'CAT',
              'Q':'CAG',
              'N':'AAC',
              'K':'AAA',
              'D':'GAT',
              'E':'GAA',
              'C':'TGC',
              'W':'TGG',
              'R':'CGT',
              'G':'GGC',
              '*':'TAA'
                  }

codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
}


sorted_codon_table = {
    'A': ['GCG','GCC','GCA','GCT'],
    'C': ['TGC','TGT'],
    'E': ['GAA','GAG'],
    'D': ['GAT','GAC'],
    'G': ['GGC','GGT','GGG','GGA'],
    'F': ['TTT','TTC'],
    'I': ['ATT','ATC','ATA'],
    'H': ['CAT','CAC'],
    'K': ['AAA','AAG'],
    'L': ['CTG','TTA','TTG','CTC','CTT','CTA'],
    'M': ['ATG'],
    'N': ['AAC','AAT'],
    'Q': ['CAG','CAA'],
    'P': ['CCG','CCA','CCT','CCC'],
    'S': ['AGC','TCG','AGT','TCC','TCT','TCA'],
    'R': ['CGC','CGT','CGG','CGA','AGA','AGG'],
    'T': ['ACC','ACG','ACT','ACA'],
    'W': ['TGG'],
    'V': ['GTG','GTT','GTC','GTA'],
    'Y': ['TAT','TAC'],
    '*': ['TAA','TGA','TAG']
    }
