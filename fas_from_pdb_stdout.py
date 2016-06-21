#!/usr/bin/env python
"""Create an (appropriately numbered) FASTA file from a PDB.

Assumes there's only one chain."""

ICODES = False

dump = ["HOH","WAT"]

trans = {"ALA":"A",
         "CYS":"C",
         "ASP":"D",
         "GLU":"E",
         "PHE":"F",
         "GLY":"G",
         "HIS":"H",
         "ILE":"I",
         "LYS":"K",
         "LEU":"L",
         "MET":"M",
         "ASN":"N",
         "PRO":"P",
         "GLN":"Q",
         "ARG":"R",
         "SER":"S",
         "THR":"T",
         "VAL":"V",
         "TRP":"W",
         "TYR":"Y",
         "MSE":"M",
        }

import sys
import os.path


def main(filename):
    CHAIN = None

    seq = {} # Map of position: (map of icode: ident)
    with open(filename,'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain = line[21]
                resn = line[17:20].upper()
                if resn in dump:
                    continue
                name = trans.get(resn,"X")
                if name == "X": #Ignore ligands/unrecognized res
                    continue
                if CHAIN is None:
                    CHAIN = chain
                if CHAIN != chain:
                     raise ValueError("Cannot deal with multiple chains: in file %s, %s/%s" % ( filename, CHAIN, chain ) )
                resi = int(line[22:26])
                icode = line[26]
                pair = (icode, name)
                rv = seq.setdefault(resi,{})
                if icode in rv:
                    if rv[icode] != name:
                        raise ValueError("Conflicting types for residue %d%s in file %s; %s/%s" % (resi, icode, filename, rv[icode], name))
                else:
                    rv[icode] = name
    if seq:
        sys.stdout.write(">"+filename+"\n")
        for i in range(1,max(seq)+1): #iterate over all indicies, even omitted ones
            rv = seq.get(i,{' ':'-'})
            icodes = rv.keys()
            icodes.sort()
            if ICODES:
                for ic in icodes:
                    sys.stdout.write(rv[ic])
            else:
                if len(icodes) != 1:
                    raise ValueError("Cannot deal with insertion codes: in file %s, position %d: '%s'" % ( filename, i, '/'.join(icodes) ) )
                sys.stdout.write(rv[icodes[0]])

        sys.stdout.write('\n')

if __name__ == "__main__":
    for filename in sys.argv[1:]:
        main(filename)
