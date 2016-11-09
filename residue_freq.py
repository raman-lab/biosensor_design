#!/usr/bin/env python
"""Outputs a table of positional frequencies (actually counts) for each
position in the input FASTA files where there is variation."""

import sys

possible = "-ACDEFGHIKLMNPQRSTVWXY"

def process_seq(seq,counts):
    for i, aa in enumerate(seq):
        pos = counts.setdefault(i+1,{})
        pos[aa] = pos.get(aa,0) + 1
        

def main(filename,counts):
    curr_seq = []
    with open(filename,'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == ">" and curr_seq:
                process_seq(''.join(curr_seq),counts)
                curr_seq = []
            elif line[0] != ">":
                curr_seq.append(line)
                
        if curr_seq:
            process_seq(''.join(curr_seq),counts)

def output(counts):
    sys.stdout.write("Pos\t" + '\t'.join(possible) + '\n')
    poses = counts.keys()
    poses.sort()
    for pos in poses:
        if len(counts[pos]) <= 1:
            continue
        sys.stdout.write(str(pos))
        for aa in possible:
            sys.stdout.write("\t" + str(counts[pos].get(aa,0)))
        sys.stdout.write('\n')
            
if __name__ == "__main__":
    counts = {}
    for filename in sys.argv[1:]:
        main(filename,counts)
    output(counts)
