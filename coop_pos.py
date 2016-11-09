#!/usr/bin/env python

"""coop_pos.py -- given a FASTA file of (aligned) top-ranked sequences,
find pairs of residues which are over represented as pairs as
compared to their individual frequencies.

Requires SciPy (for Fisher's Exact test)
"""

import sys, os
import scipy.stats

def loadseqs(fastaname):
    "Load sequences from FASTA"
    seqs = []
    with open(fastaname) as f:
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    seqs.append(''.join(seq))
                seq = []
            else:
                seq.append(line)
    if seq:
        seqs.append(''.join(seq))
    return seqs

def posfreq(seqs):
    "Computes the positional frequencies (actually counts) for all sequences."
    counts = [] # zero-based
    for seq in seqs:
        for i, aa in enumerate(seq):
            if len(counts) < i+1:
                counts.append({})
            counts[i][aa] = counts[i].get(aa,0) + 1
    return counts

def finddegenerate(counts):
    return [i for i, c in enumerate(counts) if (len(c) > 1)]

def countpairs(seqs,pos1,pos2):
    pairs = {}
    for seq in seqs:
        key = seq[pos1] + seq[pos2]
        pairs[key] = pairs.get(key,0) + 1
    return pairs
              
def main(fastaname,conf,graph=True):
    seqs = loadseqs(fastaname)
    if not seqs:
        raise ValueError("Cannot process cooperativity - no input sequences.")
    freq = posfreq(seqs)
    degen = finddegenerate(freq)
    if len(degen) < 2:
        raise ValueError("Cannot process cooperativity - need at least two positions with multiple identities.")
    degen.sort()

    graphpairs = []

    if conf is None:
        ntrials=0
        for jj, pos2 in enumerate(degen):
            for ii, pos1 in enumerate(degen):
                if ii >= jj:
                    break
                ntrials += len(freq[pos1]) * len(freq[pos2])
        conf = 0.05/ntrials; # Use Bonferroni
    
    ntrials =0
    print "pos1\tAA1\tpos2\tAA2\tobserved\texpected\tp_val"
    for jj, pos2 in enumerate(degen):
        pos2_tot = sum(freq[pos2].values())
        for ii, pos1 in enumerate(degen):
            if ii >= jj:
                break
            assert( pos1 < pos2 )
            pos1_tot = sum(freq[pos1].values())
            if pos1_tot != pos2_tot:
                raise ValueError("Sequences must be same length at varing positions.")
            pairs = countpairs(seqs,pos1,pos2)
            
            for aa1, c1 in freq[pos1].items():
                for aa2, c2 in freq[pos2].items():
                    ntrials += 1
                    key = aa1 + aa2
                    predict = pos1_tot * (float(c1)/pos1_tot) * (float(c2)/pos2_tot)
                    actual = pairs.get(key,0)
                    if (actual - 1 <= predict <= actual + 1) or (actual <= predict): #quick out for close items and negative coop
                        continue
                    table = [[actual, c1-actual], [c2-actual, pos1_tot - actual - (c1-actual) - (c2 -actual)]]
                    assert( sum([sum(i) for i in table]) == pos1_tot )
                    odds, p_val = scipy.stats.fisher_exact(table)
                    if p_val < conf:
                        graphpairs.append( ("%d%s" % (pos1 + 1, aa1), "%d%s" % (pos2 + 1, aa2), p_val ) )
                        print "%d\t%s\t%d\t%s\t%d\t%.1f\t%g" % ( pos1 + 1, aa1, pos2 + 1, aa2, actual, predict, p_val)
    sys.stderr.write("Number of trials: %d; Conf used: %0.1e; For 5%% under Bonferroni: %0.1e Sidak: %0.1e\n" % (ntrials, conf, 0.05/ntrials, 1 - (0.95)**(1.0/ntrials)) )
    if graph:
        with open("coop_graph.dot","w") as f:
            f.write("digraph G {\n")
            for pair in graphpairs:
                f.write('\t"%s" -> "%s" [dir=none];\n' % pair[:2]) # removed '"label="%0.1e"'
            f.write("}\n")
        os.system('dot -Tpdf coop_graph.dot -o coop_graph.pdf')
#os.remove("coop_graph.dot")

                               
if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "-h":
        print "Usage: coop_pos.py <fasta file name> [confidence]"
        sys.exit()
    if len(sys.argv) == 2:
        conf = None # Auto detect
    else:
        conf = float(sys.argv[2])
    main(sys.argv[1],conf)
