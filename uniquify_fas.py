#!/usr/bin/env python
"""Given one or more input FASTAs, and optionally a file of mutations to impose,
create a FASTA with only the unique sequences. Earlier FASTAs/sequences will take precident over later,
except that a later sequence with fewer mutations will take precident over an earlier one with more.

The mutation file consists of lines of mutations to impose. The line as a whole will be
imposed (if starting identities for any of the positions don't match, none of the mutations will be made.
The degeneracy symbol of 'X' is recognized as standing for any starting residue identity (e.g. X145S - 'X'
is invalid for the "to" position). If you just want to check a residue identity (and not mutate it),
omit the "to" identity (e.g. S145 ). If any of the to mutations on a line are '^', the sequence will be discarded if all the
starting residue identities on the line match.
"""

import os, sys
from optparse import OptionParser

class Mutations:
    def __init__(self,filename=None,firstpos=1):
        if filename is not None:
            self.readfromfile(filename)
        self.firstpos = firstpos

    def readfromfile(self, filename):
        self.muts = []
        with open(filename) as f:
            for line in f:
                mutset=[]
                line = line.split()
                for item in line:
                    start = item[0].upper()
                    to = item[-1].upper()
                    if start == to:
                        sys.stderr.write("WARNING: Null mutation called for with %s"%item)
                    if to == 'X':
                        sys.stderr.write("ERROR: Cannot have X as a 'to' identity.")
                        continue
                    if start == '^':
                        sys.stderr.write("ERROR: Cannot have ^ as a starting identity.")
                        continue
                    if to.isdigit():
                        pos = int(item[1:])
                        to = start
                    else:
                        pos = int(item[1:-1])
                    mutset.append( (start, pos, to) )
                self.muts.append(mutset)

    def treatseq(self, seq):
        seq = list(seq)
        nmut = 0
        for mutset in self.muts:
            for start, pos, to in mutset:
                if start != 'X' and seq[pos-self.firstpos].upper() != start:
                    break
            else: # We've made it through all of the items in the mutset, so we can make the mutations
                for start, pos, to in mutset:
                    if to == '^':  # Drop this sequence
                        return None, len(seq)
                    if seq[pos-self.firstpos].upper() != to:
                        seq[pos-self.firstpos] = to
                        nmut += 1
        return ''.join(seq), nmut

def process_seq(unique, tag, seq, muts, options):
    nmuts=0
    if muts is not None:
        seq, nmuts = muts.treatseq(seq)
        if seq is None:
            return

    if options.noblanks:
        seq = seq.replace('-','')

    options.currentpos += 1
    if seq not in unique:
        unique[seq] = (tag, nmuts, options.currentpos)
    else:
        oldtag, oldnmuts, oldpos = unique[seq]
        if oldnmuts > nmuts:
            unique[seq] = (tag, nmuts, options.currentpos)

def process_file(unique, f, options):
    tag = None
    seq = []
    for line in f:
        line = line.strip()
        if not line: #ignore blank lines
            continue
        if line.startswith(">"):
            if tag is not None:
                process_seq(unique,tag,''.join(seq),muts,options)
            seq=[]
            tag=line[1:]
        elif tag is not None:
            seq.append(line)

    if tag is not None:
        process_seq(unique,tag,''.join(seq),muts,options)

def main(args, muts, options):
    unique = {} # dictionary of post-mutation sequence:(tag, nmuts, pos)
    for name in args:
        if name == '-':
            process_file(unique, sys.stdin, options)
        else:
            with open(name) as f:
                process_file(unique, f, options)

    outseqs = []
    for seq in unique:
        tag, nmuts, pos = unique[seq]
        if options.keeporder:
            outseqs.append( (pos, tag, seq, nmuts) )
        else:
            outseqs.append( (tag, tag, seq, nmuts) )
    outseqs.sort()
    for order, tag, seq, nmuts in outseqs:
        if options.justmuts and nmuts == 0:
            continue
        if options.nomuts and nmuts != 0:
            continue
        sys.stdout.write(">%s"%tag)
        if nmuts > 0 and not options.notag:
            sys.stdout.write(" (%d mut)"%nmuts)
        sys.stdout.write("\n") # end tag line
        sys.stdout.write("%s\n"%seq)


if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [options] fasta_files ...")
    parser.add_option("-m", "--mutfile",
                  help="File listing mutations to impose")
    parser.add_option("-o", "--offset", default=1, type="int",
                  help="The first amino acid in the FASTA is position number OFFSET for the mutational identities.")
    parser.add_option("-t", "--notag", action="store_true", default=False,
                  help="Don't tag output with the number of mutations")
    parser.add_option("-n", "--nomuts", action="store_true", default=False,
                  help="Don't output sequences with mutations")
    parser.add_option("-j", "--justmuts", action="store_true", default=False,
                  help="Output just the mutated sequences (the sequences that wouldn't be outputed with -n)")
    parser.add_option("-b", "--noblanks", action="store_true", default=False,
                  help="Strip blanks ('-') from the sequences when uniquifying.")
    parser.add_option("-k", "--keeporder", action="store_true", default=False,
                  help="Keep sequences in input order.")
    options, args = parser.parse_args(sys.argv[1:])
    if not args:
        args = ['-']
    if options.mutfile is None:
        muts = None
    else:
        muts = Mutations(options.mutfile, options.offset)
    options.currentpos = 0
    try:
        main(args, muts, options)
    except IOError:
        pass

#### Old Code ########
##def parsemuts(mutfile):
##    muts = {}
##    with open(mutfile) as f:
##        for line in f:
##            line = line.split()
##            for item in line:
##                start = item[0].upper()
##                to = item[-1].upper()
##                if start == to:
##                    sys.stderr.write("WARNING: Null mutation called for with %s"%item)
##                    continue
##                pos = int(item[1:-1])
##                pos_dict = muts.setdefault(pos,{})
##                if start in pos_dict and pos_dict[start] != to:
##                    raise ValueError("Mutation file gives inconsistent mutations for %s%d: %s and %s"%(start,pos,pos_dict[start],to))
##                pos_dict[start] = to
##    return muts
##
##def process_seq(unique,tag,seq,muts):
##    nmuts=0
##    if muts:
##        seq=list(seq)
##        for pos in muts:
##            aa = seq[pos-1] #zero vs 1 based
##            if aa in muts[pos]:
##                seq[pos-1] = muts[pos][aa] #zero vs 1 based
##                nmuts += 1
##            elif 'X' in muts[pos] and aa != muts[pos]['X']:
##                seq[pos-1] = muts[pos]['X'] #zero vs 1 based
##                nmuts += 1
##        seq=''.join(seq)
##    if seq not in unique:
##        unique[seq] = (tag, nmuts)
##    else:
##        oldtag, oldnmuts = unique[seq]
##        if nmuts < oldnmuts:
##            unique[seq] = (tag, nmuts)


