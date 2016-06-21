#! /usr/bin/python

from phil import *
from Bio import SeqIO
import stats


def Help():
    print "Computes the number the mutations that in Rosetta design sequence"
    print "Usage: <lacI_WT.fasta> <library of sequences> <start pos> <end pos>"
    exit()

if len(argv) < 2:
    Help()


lacI_WT = sys.argv[1]
design_sequences = sys.argv[2]
start_pos = int(sys.argv[3])
end_pos = int(sys.argv[4])


#for record1 in SeqIO.parse(handle1, "fasta"):
#    print record1.seq
#    for i in range(start_pos-1, end_pos):
#        print record1.seq[i]
#       handle2 = open(design_sequences, "rU")
#        for record2 in SeqIO.parse(handle2,"fasta"):
#            for j in range(0,len(record2.seq)):
#                print record1.seq[i], record2.seq[j]

mutation_freq_per_pos = {}

for i in range(start_pos,end_pos+1):
    mutation_freq_per_pos[i] = 0
#    print i

#print mutation_freq_per_pos

def hamdist(WT,query):
    u = zip(WT, query)
    d = dict(u)
    y = []
    counter = 0
    pos = 0
    for i,j in u:
        actual_pos = start_pos + pos
        pos += 1
#        print "i, j ", i, j
        if i == j:
            continue
        if i!=j:
#            print i,actual_pos,j
            sys.stdout.write('%s%s%s\n'%(i,actual_pos,j))
            counter += 1
#            print "actual pos", actual_pos
            mutation_freq_per_pos[actual_pos] += 1

    return counter

num_diff = []
handle1 = open(lacI_WT, "rU")
total_sequences = 0
for record1 in SeqIO.parse(handle1,"fasta"):
    WT = record1.seq[start_pos-1:end_pos]
    handle2 = open(design_sequences,"rU")
    for record2 in SeqIO.parse(handle2,"fasta"):
#        print WT, record2.seq
        num_diff.append(hamdist(WT, record2.seq))
        total_sequences += 1
        print "**********"
 
exit()
#print num_diff

numbins = 10
lowest_bin = 0
highest_bin = 10
mutation_dist = stats.lhistogram(num_diff,numbins,(lowest_bin,highest_bin),0)[0]
bin_width = stats.lhistogram(num_diff,numbins,(lowest_bin,highest_bin),0)[2]

print stats.lhistogram(num_diff,numbins,(lowest_bin,highest_bin),0)

#print mutation_dist
#print bin_width
for i in range(lowest_bin, highest_bin):
    sys.stdout.write('%s\t%s\n'%(i, mutation_dist[i]))

#print mutation_freq_per_pos

for keys in mutation_freq_per_pos:
    WT_res = record1.seq[keys-1]
    if mutation_freq_per_pos[keys] != 0:
#        print keys, mutation_freq_per_pos[keys], record1.seq[keys-1]
        sys.stdout.write('%s%s %2.2f\n'%(WT_res,keys, float(mutation_freq_per_pos[keys])/total_sequences))

print "total seq ", total_sequences
