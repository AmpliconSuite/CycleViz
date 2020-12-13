#!/usr/bin/env python

import argparse
import os
import subprocess

try:
    from subprocess import DEVNULL  # Python 3.
except ImportError:
    DEVNULL = open(os.devnull, 'wb')


# run a command and yield the stdout on the fly, derived from user "tokland" on StackOverflow
# https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
def execute(cmd, redirect_stderr = False):
    stderr_loc = None
    if redirect_stderr:
        stderr_loc = DEVNULL

    popen = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True, stderr=stderr_loc)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


# MAIN

#inputs: bam file, bed file
parser = argparse.ArgumentParser(description="Get coverage bedgraph and baseline coverage")
parser.add_argument("--bam",type=str, help="path to bam file", required=True)
parser.add_argument("--bed",type=str, help="path to bed file", required=True)
parser.add_argument("--estimate_average_coverage", help="get an estimate of the average coverage per chromosome",
                    action='store_true')
parser.add_argument("-o", type=str,help="output file prefix. Default to prefix of input bed")

args = parser.parse_args()

with open(args.bed) as infile:
    regions = []
    for line in infile:
        fields = line.rstrip().rsplit()
        regions.append(fields)

base_cmd = "samtools mpileup -B -d 50000 -q 5 -r "
bambase = os.path.splitext(os.path.basename(args.bam))[0]
if not args.o:
    bedbase = os.path.splitext(os.path.basename(args.bed))[0]
    bgout = bedbase + "_position_coverage.bedgraph"

else:
    bgout = args.o + "_position_coverage.bedgraph"

# make the bedgraph
with open(bgout, 'w') as outfile:
    for ind, r in enumerate(regions):
        prevpoint = ["",0,0,0]
        cmd = base_cmd + r[0] + ":" + r[1] + "-" + r[2] + " " + args.bam + " | cut -f 1-2,4"
        for pline in execute(cmd):
            fields = pline.rstrip().rsplit()
            c, p, n = fields[0], int(fields[1]), int(fields[2])
            if c != prevpoint[0] or n != prevpoint[3] or p > prevpoint[2] + 1:
                if prevpoint[0]:
                    outfile.write("\t".join([str(x) for x in prevpoint]) + "\n")

                prevpoint = [c, p, p, n]

            else:
                prevpoint[2] = p

        if prevpoint[0]:
            outfile.write("\t".join([str(x) for x in prevpoint]) + "\n")


# get a bedgraph of the mean coverage values
if args.estimate_average_coverage:
    # first get the average read length
    # command based on solution by Chris Miller at BioStars. https://www.biostars.org/p/65216/
    print("estimating read length")
    cmd = "samtools view -F 4 " + args.bam + " | head -n 2000000 | cut -f 10 | perl -ne \'chomp;print length($_) . " + "\"" + "\\n" + "\"" + "' | sort | uniq -c"
    print(cmd)
    runT = 0.0
    runW = 0.0
    for i in execute(cmd, redirect_stderr=True):
        c, l = [int(x) for x in i.rstrip().lstrip().rsplit()]
        runT+=c
        runW+=(c*l)

    meanRL = runW/runT

    # now get the number of aligned reads per chromosome, the length of the chr and do Lander-Waterman stats
    print("getting chromosome coverage")
    cmd = "samtools idxstats " + args.bam
    print(cmd)
    with open(bambase + "_chromosome_coverage.bedgraph",'w') as outfile:
        for covline in execute(cmd):
            fields = covline.rstrip().rsplit()
            if fields[0] != '*':
                G = int(fields[1])
                n = int(fields[2])
                a = 0
                if G > 0:
                    a = n*meanRL/G

                outfile.write("\t".join([fields[0], "1", fields[1], str(a)]) + "\n")
