#!/usr/bin/env python

import os
import sys

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

rdd = {'+':'-', '-':'+'}

# output a reversed alignment from a segaligner _aln.txt file format.

alignment_file = sys.argv[1]

oname = os.path.splitext(os.path.basename(alignment_file))[0]
oname = oname.rstrip("\_aln")
oname += "_reversed_aln.txt"

with open(alignment_file) as infile, open(oname, 'w') as ofile:
    ofile.write(next(infile))
    seqline = next(infile)

    # reverse the seq
    fields = seqline.rsplit("\t")
    seq = fields[0][1:].rsplit(",")
    rseq = [x[:-1] + "-" if x[-1] == "+" else x[:-1] + "+" for x in seq[::-1]]
    fields[0] = "#" + ",".join(rseq)
    ofile.write("\t".join(fields))

    # header
    ofile.write(next(infile))

    # alignment lines
    lines = [l for l in infile]
    max_ind = len(seq) - 1
    if seq[0] == "0+":
        leaveInitInd = False
        max_ind -= 2

    else:
        leaveInitInd = True

    rev_lines = []
    for x in lines[::-1]:
        fields = x.rsplit("\t")
        fields[4] = rdd[fields[4]]
        fields[5] = rdd[fields[5]]
        san = int(fields[6])
        if san != 0 or not leaveInitInd:
            fields[6] = str(max_ind - san)

        rev_lines.append("\t".join(fields))

    for oline in rev_lines:
        ofile.write(oline)
