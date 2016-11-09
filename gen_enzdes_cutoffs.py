#!/usr/bin/env python
"""Generates an enzdes cutoff file based on the current tabular values.

Currently uses median absolute deviation to determine cutoffs.

Can take multiple input tables, and requires a descriptor file:
Use "%(LIG)s" to refer to the ligans descriptor (e.g. "SR_3")
Line is "value column" "keep less/greater than" "cutoff value above (neg:below) center" "Type of descriptor"
MAD - median absolute deviation from median
% - percentile (positive starts from lower end (0=min), negative starts from upper end (0=max))
SD - standard deviation from mean
ABS - absolute value

---Example descriptor---
total_score < 2 MAD
fa_rep < 2 MAD

%(LIG)s_total_score < 2 MAD
%(LIG)s_interf_E_1_2 < -1 MAD
------------


MAD table:
----------
~1.4826 times the MAD is the SD for normally distributed data
(Both are linearly related as they are in the measurement units).
The fraction within z in either direction of the mean is erf( z/sqrt(2) )
so the one tailed outliers are (1 - erf( z/sqrt(2) ))/2

MAD     SD eq (normal)      One tailed % (normal)
1       0.67                25
1.48    1                   15.9
2       1.35                 8.9
2.5     1.69                 4.6
3       2.02                 2.2

"""

import os, sys
from optparse import OptionParser
from math import sqrt

def parse_descriptors(descriptorfile):
    tags = []
    filters = []
    with open(descriptorfile) as f:
        for line in f:
            line = line.split()
            if (not line) or (line[0].startswith("#")):
                continue
            tags.append(line[0])
            filters.append(line[0:4])
    return tags, filters

def parse_tables(inputfiles, tags):
    entries = []
    output_subs = None
    for name in inputfiles:
        new_tags = []
        with open(name) as f:
            header = f.readline().split()
            SRs = [h for h in header if h.startswith("SR")]
            if SRs:
                subs = {"LIG":SRs[-1][0:4]} # Doesn't work with more than 9 constraints, but I'm punting
            else:
                subs = {}
            if output_subs is None:
                output_subs = subs
            if output_subs != subs:
                raise ValueError("Inconsistent substitution patterns between multiple files: %s, %s"%(output_subs, subs))
            for t in tags:
                if t%subs not in header:
                    raise ValueError("Could not find tag '%s' in header of '%s': %s" %(t%subs, name, ' '.join(header))) 
            indexes = [header.index(t%subs) for t in tags]
            for line in f:
                line = line.split()
                if (not line) or (line[0].startswith("#")) or (line[0][0].isalpha()):
                    continue
                entries.append([line[i] for i in indexes])
    return output_subs, entries

def ABS_cutoff(items, n, op, dev):
    return float(dev)

def MEDIAN(seq):
    values = list(seq)
    values.sort()
    nvals = len(values)
    if nvals%2 == 1:
        median = values[ int(nvals/2) ]
    else:
        median = (values[ int(nvals/2)-1 ] + values[ int(nvals/2) ])/2
    return median
    
def MAD_cutoff(items, n, op, dev):
    values = [float(i[n]) for i in items]
    median = MEDIAN(values)
    abs_dev = [ abs( v - median ) for v in values ]
    mad = MEDIAN(abs_dev)
    cutoff = median + (float(dev) * mad)
    #Try to account for rounding issues
    if op == "<":
        cutoff += 0.01
    if op == ">":
        cutoff -= 0.01
    return cutoff

def percentile_cutoff(items, n, op, dev):
    values = [float(i[n]) for i in items]
    values.sort()
    pos = float(dev)/100*len(values) - 0.5 # 0.5 is to adjust for zero indexing issues
    if pos == int(pos):
        cutoff = values[int(pos)]
    elif pos < 0:
        cutoff = values[0]
    elif pos > len(values) - 1:
        cutoff = values[-1]
    else:
        frac = pos - int(pos)
        cutoff = (1-frac)*values[int(pos)] + (frac)*values[int(pos)+1] #Fancy linear interpolation
    #Try to account for rounding issues
    if op == "<":
        cutoff += 0.01
    if op == ">":
        cutoff -= 0.01
    return cutoff
        

def SD_cutoff(items, n, op, dev):
    values = [float(i[n]) for i in items]
    mean = sum(values)/len(values)
    sd = sqrt( sum([ (v-mean)**2 for v in values ])/len(values) )
    cutoff = mean + (float(dev) * sd)
    #Try to account for rounding issues
    if op == "<":
        cutoff += 0.01
    if op == ">":
        cutoff -= 0.01
    return cutoff
    
def find_cutoffs(items, filters):
    cutoff_vals = []
    for n, filter in enumerate(filters):
        tag, op, dev, type = filter
        type = type.upper()
        if type=="ABS":
            cutoff_vals.append(ABS_cutoff(items, n, op, dev))
        elif type=="MAD":
            cutoff_vals.append(MAD_cutoff(items, n, op, dev))
        elif type=="%":
            cutoff_vals.append(percentile_cutoff(items, n, op, dev))
        elif type=="SD":
            cutoff_vals.append(SD_cutoff(items, n, op, dev))
        else:
            raise ValueError("Descriptor type '%s' not recognized for line: '%s'"%(type,' '.join(filters)))
    assert(len(cutoff_vals) == len(filters))
    return cutoff_vals

def write_cutoff_file(outfile, filters, cutoffs, subs):
    with open(outfile,'w') as f:
        for n, filter in enumerate(filters):
            tag, op, dev, type = filter
            f.write("req %s value %s %.2f\n"%(tag%subs, op, cutoffs[n]))

def main(inputfiles, descriptorfile, outfile):
    tags, filters = parse_descriptors(descriptorfile)
    subs, items = parse_tables(inputfiles, tags)
    cutoffs = find_cutoffs(items, filters)
    write_cutoff_file(outfile, filters, cutoffs, subs) 
    
if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [options] enzdes_table_files ...")
    parser.add_option("-c", "--cutoff",
                  help="Descriptor specification for cutoff")
    parser.add_option("-o", "--out", default="cutoffs",
                  help="Filename to output to")
    options, args = parser.parse_args(sys.argv[1:])
    test = parser.parse_args()
    if not args:
        raise ValueError("Must specify at least one input tables")
    if options.cutoff is None:
        raise ValueError("Must specify descriptor specification")
    main(args, options.cutoff, options.out)
        
    
