#!/usr/bin/env python

# Reverse the path/cycles in an AA cycles file and write a new file named *_reversed_cycles.txt.
# E.g. turn 4+,5+,2- into 2+,5-,4-. Reversal applied to each entry in the cycles file.
# Optionally can also reverse an AR alignment file.

"""
Jens Luebeck
UC San Diego, Bioinformatics & Systems Biology
jluebeck@ucsd.edu
"""

import argparse
import os
import subprocess

CV_RESOURCES = os.path.dirname(os.path.abspath(__file__))
rdd = {'+':'-', '-':'+'}


def reverse_AA_cycles_file(cf, ofname):
    with open(cf) as infile, open(ofname, 'w') as outfile:
        for line in infile:
            if not line.startswith("Cycle"):
                outfile.write(line)
            else:
                cfields = line.rstrip().rsplit(";")
                cvals = [(x.rsplit("=")[0], x.rsplit("=")[1]) for x in cfields]
                for ind, (x, y) in enumerate(cvals):
                    if x.startswith("Segments"):
                        sstring = y.rsplit(",")
                        slist = [(s[:-1], s[-1]) for s in sstring]
                        revslist = [s[0] + rdd[s[1]] for s in slist[::-1]]
                        revy = ",".join(revslist)
                        cvals[ind] = (x, revy)

                rcfields = [x + "=" + y for x, y in cvals]
                outfile.write(";".join(rcfields) + "\n")

# MAIN
# inputs: bam file, bed file
parser = argparse.ArgumentParser(description="Reverse one or more entries in a cycles file.")
parser.add_argument("-c", type=str, help="path to cycles file", required=True)
parser.add_argument("--aln_file", type=str, help="path to AR OM alignment file (*_aln.txt) (optional)")

args = parser.parse_args()

ofname = args.c.rstrip("cycles.txt") + "reversed_cycles.txt"
reverse_AA_cycles_file(args.c, ofname)
print("Reversed AA cycles file " + args.c)


if args.aln_file:
    cmd = CV_RESOURCES + "/reverse_aln_file.py " + args.aln_file
    print(cmd)
    subprocess.call(cmd, shell=True)
    print("Reversed AR alignment file")

